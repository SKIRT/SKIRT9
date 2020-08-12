/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "MediumSystem.hpp"
#include "Configuration.hpp"
#include "Constants.hpp"
#include "DensityInCellInterface.hpp"
#include "DisjointWavelengthGrid.hpp"
#include "FatalError.hpp"
#include "LockFree.hpp"
#include "Log.hpp"
#include "LyaUtils.hpp"
#include "MaterialMix.hpp"
#include "MediumState.hpp"
#include "NR.hpp"
#include "Parallel.hpp"
#include "ParallelFactory.hpp"
#include "PathSegmentGenerator.hpp"
#include "PhotonPacket.hpp"
#include "ProcessManager.hpp"
#include "Random.hpp"
#include "ShortArray.hpp"
#include "StringUtils.hpp"

////////////////////////////////////////////////////////////////////

namespace
{
    // maximum number of cell densities calculated between two invocations of infoIfElapsed()
    const size_t logProgressChunkSize = 10000;
}

////////////////////////////////////////////////////////////////////

void MediumSystem::setupSelfAfter()
{
    SimulationItem::setupSelfAfter();
    auto log = find<Log>();
    auto parfac = find<ParallelFactory>();
    _config = find<Configuration>();

    // ----- allocate memory -----

    _numCells = _grid->numCells();
    if (_numCells < 1) throw FATALERROR("The spatial grid must have at least one cell");
    _numMedia = _media.size();

    // initial state
    size_t allocatedBytes = 0;
    _state1v.resize(_numCells);
    allocatedBytes += _state1v.size() * sizeof(State1);
    _state2vv.resize(_numCells * _numMedia);
    allocatedBytes += _state2vv.size() * sizeof(State2);

    // radiation field
    if (_config->hasRadiationField())
    {
        _wavelengthGrid = _config->radiationFieldWLG();
        _rf1.resize(_numCells, _wavelengthGrid->numBins());
        allocatedBytes += _rf1.size() * sizeof(double);

        if (_config->hasSecondaryRadiationField())
        {
            _rf2.resize(_numCells, _wavelengthGrid->numBins());
            _rf2c.resize(_numCells, _wavelengthGrid->numBins());
            allocatedBytes += 2 * _rf2.size() * sizeof(double);
        }
    }

    // inform user
    log->info(typeAndName() + " allocated " + StringUtils::toMemSizeString(allocatedBytes) + " of memory");

    // ----- calculate cell densities, bulk velocities, and volumes in parallel -----

    log->info("Calculating densities for " + std::to_string(_numCells) + " cells...");
    auto dic = _grid->interface<DensityInCellInterface>(0, false);  // optional fast-track interface for densities
    int numSamples = _config->numDensitySamples();
    bool oligo = _config->oligochromatic();
    int hMag = _config->magneticFieldMediumIndex();
    int hLya = _config->lyaMediumIndex();
    log->infoSetElapsed(_numCells);
    parfac->parallelDistributed()->call(
        _numCells, [this, log, dic, numSamples, oligo, hMag, hLya](size_t firstIndex, size_t numIndices) {
            ShortArray<8> nsumv(_numMedia);

            while (numIndices)
            {
                size_t currentChunkSize = min(logProgressChunkSize, numIndices);
                for (size_t m = firstIndex; m != firstIndex + currentChunkSize; ++m)
                {
                    Position center = _grid->centralPositionInCell(m);

                    // density: use optional fast-track interface or sample 100 random positions within the cell
                    if (dic)
                    {
                        for (int h = 0; h != _numMedia; ++h) state(m, h).n = dic->numberDensity(h, m);
                    }
                    else
                    {
                        nsumv.clear();
                        for (int n = 0; n < numSamples; n++)
                        {
                            Position bfr = _grid->randomPositionInCell(m);
                            for (int h = 0; h != _numMedia; ++h) nsumv[h] += _media[h]->numberDensity(bfr);
                        }
                        for (int h = 0; h != _numMedia; ++h) state(m, h).n = nsumv[h] / numSamples;
                    }

                    // bulk velocity: weighted average at cell center; for oligochromatic simulations, leave at zero
                    if (!oligo)
                    {
                        double n = 0.;
                        Vec v;
                        for (int h = 0; h != _numMedia; ++h)
                        {
                            n += state(m, h).n;
                            v += state(m, h).n * _media[h]->bulkVelocity(center);
                        }
                        if (n > 0.) state(m).v = v / n;  // leave bulk velocity at zero if cell has no material
                    }

                    // magnetic field: retrieve from medium component that specifies it, if any
                    if (hMag >= 0)
                    {
                        state(m).B = _media[hMag]->magneticField(center);
                    }

                    // gas temperature: retrieve from medium component that specifies it, if any
                    if (hLya >= 0)
                    {
                        // leave the temperature at zero if the cell does not contain any Lya gas;
                        // otherwise make sure the temperature is at least the local universe CMB temperature
                        if (state(m, hLya).n > 0.)
                            state(m).T = max(Constants::Tcmb(), _media[hLya]->temperature(center));
                    }

                    // volume
                    state(m).V = _grid->volume(m);
                }
                log->infoIfElapsed("Calculated cell densities: ", currentChunkSize);
                firstIndex += currentChunkSize;
                numIndices -= currentChunkSize;
            }
        });

    // communicate the calculated states across multiple processes, if needed
    communicateStates();

    log->info("Done calculating cell densities");

    // ----- obtain the material mix pointers -----

    for (int m = 0; m != _numCells; ++m)
    {
        Position bfr = _grid->centralPositionInCell(m);
        for (int h = 0; h != _numMedia; ++h) state(m, h).mix = _media[h]->mix(bfr);
    }
}

////////////////////////////////////////////////////////////////////

void MediumSystem::communicateStates()
{
    if (!ProcessManager::isMultiProc()) return;

    // NOTE: once the design of the state data structures is stable, a custom communication procedure could be provided
    //       in the meantime, we copy the data into a temporary table so we can use the standard sumToAll procedure
    Table<2> data;

    // volumes, bulk velocities, and magnetic fields
    data.resize(_numCells, 8);
    for (int m = 0; m != _numCells; ++m)
    {
        data(m, 0) = state(m).V;
        data(m, 1) = state(m).v.x();
        data(m, 2) = state(m).v.y();
        data(m, 3) = state(m).v.z();
        data(m, 4) = state(m).B.x();
        data(m, 5) = state(m).B.y();
        data(m, 6) = state(m).B.z();
        data(m, 7) = state(m).T;
    }
    ProcessManager::sumToAll(data.data());
    for (int m = 0; m != _numCells; ++m)
    {
        state(m).V = data(m, 0);
        state(m).v = Vec(data(m, 1), data(m, 2), data(m, 3));
        state(m).B = Vec(data(m, 4), data(m, 5), data(m, 6));
        state(m).T = data(m, 7);
    }

    // densities
    data.resize(_numCells, _numMedia);
    for (int m = 0; m != _numCells; ++m)
        for (int h = 0; h != _numMedia; ++h) data(m, h) = state(m, h).n;
    ProcessManager::sumToAll(data.data());
    for (int m = 0; m != _numCells; ++m)
        for (int h = 0; h != _numMedia; ++h) state(m, h).n = data(m, h);
}

////////////////////////////////////////////////////////////////////

int MediumSystem::dimension() const
{
    int result = 1;
    for (auto medium : _media) result = max(result, medium->dimension());
    return result;
}

////////////////////////////////////////////////////////////////////

int MediumSystem::gridDimension() const
{
    return _grid->dimension();
}

////////////////////////////////////////////////////////////////////

int MediumSystem::numMedia() const
{
    return _numMedia;
}

////////////////////////////////////////////////////////////////////

int MediumSystem::numCells() const
{
    return _numCells;
}

////////////////////////////////////////////////////////////////////

double MediumSystem::volume(int m) const
{
    return state(m).V;
}

////////////////////////////////////////////////////////////////////

Vec MediumSystem::bulkVelocity(int m) const
{
    return state(m).v;
}

////////////////////////////////////////////////////////////////////

Vec MediumSystem::magneticField(int m) const
{
    return state(m).B;
}

////////////////////////////////////////////////////////////////////

double MediumSystem::gasTemperature(int m) const
{
    return state(m).T;
}

////////////////////////////////////////////////////////////////////

bool MediumSystem::hasMaterialType(MaterialMix::MaterialType type) const
{
    for (int h = 0; h != _numMedia; ++h)
        if (state(0, h).mix->materialType() == type) return true;
    return false;
}

////////////////////////////////////////////////////////////////////

bool MediumSystem::isMaterialType(MaterialMix::MaterialType type, int h) const
{
    return state(0, h).mix->materialType() == type;
}

////////////////////////////////////////////////////////////////////

double MediumSystem::numberDensity(int m, int h) const
{
    return state(m, h).n;
}

////////////////////////////////////////////////////////////////////

double MediumSystem::massDensity(int m, int h) const
{
    return state(m, h).n * state(m, h).mix->mass();
}

////////////////////////////////////////////////////////////////////

const MaterialMix* MediumSystem::mix(int m, int h) const
{
    return state(m, h).mix;
}

////////////////////////////////////////////////////////////////////

double MediumSystem::opacityAbs(double lambda, int m, int h) const
{
    MediumState mst(this, m, h);
    return state(m, h).mix->opacityAbs(lambda, &mst, nullptr);
}

////////////////////////////////////////////////////////////////////

double MediumSystem::opacitySca(double lambda, int m, int h) const
{
    MediumState mst(this, m, h);
    return state(m, h).mix->opacitySca(lambda, &mst, nullptr);
}

////////////////////////////////////////////////////////////////////

double MediumSystem::opacityExt(double lambda, int m, int h) const
{
    MediumState mst(this, m, h);
    return state(m, h).mix->opacityExt(lambda, &mst, nullptr);
}

////////////////////////////////////////////////////////////////////

double MediumSystem::opacityAbs(double lambda, int m, MaterialMix::MaterialType type) const
{
    double result = 0.;
    for (int h = 0; h != _numMedia; ++h)
        if (state(0, h).mix->materialType() == type) result += opacityAbs(lambda, m, h);
    return result;
}

////////////////////////////////////////////////////////////////////

double MediumSystem::opacityExt(double lambda, int m, MaterialMix::MaterialType type) const
{
    double result = 0.;
    for (int h = 0; h != _numMedia; ++h)
        if (state(0, h).mix->materialType() == type) result += opacityExt(lambda, m, h);
    return result;
}

////////////////////////////////////////////////////////////////////

double MediumSystem::opacityExt(double lambda, int m) const
{
    double result = 0.;
    for (int h = 0; h != _numMedia; ++h) result += opacityExt(lambda, m, h);
    return result;
}

////////////////////////////////////////////////////////////////////

double MediumSystem::albedo(double lambda, int m) const
{
    double ksca = 0.;
    double kext = 0.;
    for (int h = 0; h != _numMedia; ++h)
    {
        ksca += opacitySca(lambda, m, h);
        kext += opacityExt(lambda, m, h);
    }
    return kext > 0. ? ksca / kext : 0.;
}

////////////////////////////////////////////////////////////////////

double MediumSystem::perceivedWavelengthForScattering(const PhotonPacket* pp)
{
    if (_config->hasMovingMedia())
        return pp->perceivedWavelength(state(pp->interactionCellIndex()).v,
                                       _config->lyaExpansionRate() * pp->interactionDistance());
    else
        return pp->wavelength();
}

////////////////////////////////////////////////////////////////////

bool MediumSystem::weightsForScattering(Array& wv, double lambda, const PhotonPacket* pp)
{
    // resize the target array
    wv.resize(_numMedia);

    // for a single component, the weight is trivial
    if (_numMedia == 1)
    {
        wv[0] = 1.;
        return true;
    }

    // locate the cell hosting the scattering event
    int m = pp->interactionCellIndex();

    // calculate the weights and their sum
    double sum = 0.;
    for (int h = 0; h != _numMedia; ++h)
    {
        MediumState mst(this, m, h);
        wv[h] = state(m, h).mix->opacitySca(lambda, &mst, pp);
        sum += wv[h];
    }

    // normalize the weights
    if (sum > 0)
    {
        wv /= sum;
        return true;
    }

    // none of the media scatter
    return false;
}

////////////////////////////////////////////////////////////////////

void MediumSystem::peelOffScattering(double lambda, const Array& wv, Direction bfkobs, Direction bfky, PhotonPacket* pp,
                                     PhotonPacket* ppp)
{
    // get the cell hosting the scattering event
    int m = pp->interactionCellIndex();

    // calculate the weighted sum of the effects on the Stokes vector and on the wavelength for all media
    double I = 0., Q = 0., U = 0., V = 0.;
    for (int h = 0; h != _numMedia; ++h)
    {
        MediumState mst(this, m, h);
        state(m, h).mix->peeloffScattering(I, Q, U, V, lambda, wv[h], bfkobs, bfky, &mst, pp);
    }

    // pass the result to the peel-off photon packet
    ppp->launchScatteringPeelOff(pp, bfkobs, bulkVelocity(m), lambda, I);
    if (_config->hasPolarization()) ppp->setPolarized(I, Q, U, V, pp->normal());
}

////////////////////////////////////////////////////////////////////

void MediumSystem::simulateScattering(Random* random, PhotonPacket* pp)
{
    // locate the cell hosting the scattering event
    int m = pp->interactionCellIndex();

    // calculate the perceived wavelength in the cell
    double lambda = perceivedWavelengthForScattering(pp);

    // select a medium component within that cell
    int h = 0;
    if (_numMedia > 1)
    {
        // build the cumulative distribution corresponding to the scattering opacities
        Array Xv;
        NR::cdf(Xv, _numMedia, [this, lambda, pp, m](int h) {
            MediumState mst(this, m, h);
            return state(m, h).mix->opacitySca(lambda, &mst, pp);
        });

        // randomly select an index
        h = NR::locateClip(Xv, random->uniform());
    }

    // actually perform the scattering event for this cell and medium component
    MediumState mst(this, m, h);
    state(m, h).mix->performScattering(lambda, &mst, pp);
}

////////////////////////////////////////////////////////////////////

namespace
{
    // This function returns a thread-local instance of the path segment generator for the specified grid
    // that is initialized to the starting position and direction of the specified path.
    // Providing a thread-local instance avoids creating a new generator for each use.
    PathSegmentGenerator* getPathSegmentGenerator(SpatialGrid* grid, const SpatialGridPath* path)
    {
        thread_local SpatialGrid* t_grid{nullptr};
        thread_local std::unique_ptr<PathSegmentGenerator> t_generator;

        if (grid != t_grid)
        {
            t_grid = grid;
            t_generator = grid->createPathSegmentGenerator();
        }
        t_generator->start(path);
        return t_generator.get();
    }
}

////////////////////////////////////////////////////////////////////

double MediumSystem::getOpticalDepth(const SpatialGridPath* path, double lambda, MaterialMix::MaterialType type)
{
    // determine the geometric details of the path and calculate the optical depth at the same time
    auto generator = getPathSegmentGenerator(_grid, path);
    double tau = 0.;
    while (generator->next())
    {
        if (generator->m() >= 0) tau += opacityExt(lambda, generator->m(), type) * generator->ds();
    }
    return tau;
}

////////////////////////////////////////////////////////////////////

void MediumSystem::setOpticalDepths(PhotonPacket* pp)
{
    // determine and store the path segments in the photon packet
    auto generator = getPathSegmentGenerator(_grid, pp);
    pp->clear();
    while (generator->next())
    {
        pp->addSegment(generator->m(), generator->ds());
    }

    // calculate the cumulative optical depth and store it in the photon packet for each path segment
    double tau = 0.;
    int i = 0;

    // single medium, spatially constant cross sections
    if (_config->hasSingleConstantMedium())
    {
        double section = state(0, 0).mix->sectionExt(pp->wavelength());
        for (auto& segment : pp->segments())
        {
            if (segment.m >= 0) tau += section * state(segment.m, 0).n * segment.ds;
            pp->setOpticalDepth(i++, tau);
        }
    }

    // multiple media, spatially constant cross sections
    else if (_config->hasMultipleConstantMedia())
    {
        ShortArray<8> sectionv(_numMedia);
        for (int h = 0; h != _numMedia; ++h) sectionv[h] = state(0, h).mix->sectionExt(pp->wavelength());
        for (auto& segment : pp->segments())
        {
            if (segment.m >= 0)
                for (int h = 0; h != _numMedia; ++h) tau += sectionv[h] * state(segment.m, h).n * segment.ds;
            pp->setOpticalDepth(i++, tau);
        }
    }

    // spatially variable cross sections
    else
    {
        for (auto& segment : pp->segments())
        {
            if (segment.m >= 0)
            {
                double lambda = pp->perceivedWavelength(state(segment.m).v, _config->lyaExpansionRate() * segment.s);
                tau += opacityExt(lambda, segment.m) * segment.ds;
            }
            pp->setOpticalDepth(i++, tau);
        }
    }
}

////////////////////////////////////////////////////////////////////

bool MediumSystem::setInteractionPoint(PhotonPacket* pp, double tauscat)
{
    auto generator = getPathSegmentGenerator(_grid, pp);
    double tau = 0.;
    double s = 0.;

    // loop over the segments of the path until the interaction optical depth is reached or the path ends

    // --> single medium, spatially constant cross sections
    if (_config->hasSingleConstantMedium())
    {
        double section = state(0, 0).mix->sectionExt(pp->wavelength());
        while (generator->next())
        {
            // remember the cumulative optical depth and distance at the start of this segment
            // so that we can interpolate the interaction point should it happen to be inside this segment
            double tau0 = tau;
            double s0 = s;

            // calculate the cumulative optical depth and distance at the end of this segment
            double ds = generator->ds();
            int m = generator->m();
            if (m >= 0) tau += section * state(m, 0).n * generator->ds();
            s += ds;

            // if the interaction point is inside this segment, store it in the photon packet
            if (tauscat < tau)
            {
                pp->setInteractionPoint(m, NR::interpolateLinLin(tauscat, tau0, tau, s0, s));
                return true;
            }
        }
    }

    // --> multiple media, spatially constant cross sections
    else if (_config->hasMultipleConstantMedia())
    {
        ShortArray<8> sectionv(_numMedia);
        for (int h = 0; h != _numMedia; ++h) sectionv[h] = state(0, h).mix->sectionExt(pp->wavelength());
        while (generator->next())
        {
            // remember the cumulative optical depth and distance at the start of this segment
            // so that we can interpolate the interaction point should it happen to be inside this segment
            double tau0 = tau;
            double s0 = s;

            // calculate the cumulative optical depth and distance at the end of this segment
            double ds = generator->ds();
            int m = generator->m();
            if (m >= 0)
                for (int h = 0; h != _numMedia; ++h) tau += sectionv[h] * state(m, h).n * ds;
            s += ds;

            // if the interaction point is inside this segment, store it in the photon packet
            if (tauscat < tau)
            {
                pp->setInteractionPoint(m, NR::interpolateLinLin(tauscat, tau0, tau, s0, s));
                return true;
            }
        }
    }

    // --> spatially variable cross sections
    else
    {
        while (generator->next())
        {
            // remember the cumulative optical depth and distance at the start of this segment
            // so that we can interpolate the interaction point should it happen to be inside this segment
            double tau0 = tau;
            double s0 = s;

            // calculate the cumulative optical depth and distance at the end of this segment
            double ds = generator->ds();
            int m = generator->m();
            if (m >= 0)
            {
                double lambda = pp->perceivedWavelength(state(m).v, _config->lyaExpansionRate() * s);
                tau += opacityExt(lambda, m) * ds;
            }
            s += ds;

            // if the interaction point is inside this segment, store it in the photon packet
            if (tauscat < tau)
            {
                pp->setInteractionPoint(m, NR::interpolateLinLin(tauscat, tau0, tau, s0, s));
                return true;
            }
        }
    }

    // the interaction point is outside of the path
    return false;
}

////////////////////////////////////////////////////////////////////

double MediumSystem::getOpticalDepth(PhotonPacket* pp, double distance)
{
    // determine the optical depth at which the packet's contribution becomes zero
    // or abort right away if the contribution is zero to begin with
    double L = pp->luminosity();
    if (L <= 0) return std::numeric_limits<double>::infinity();
    double taumax = std::log(L) + 745;

    // determine the geometric details of the path and calculate the optical depth at the same time
    auto generator = getPathSegmentGenerator(_grid, pp);
    double tau = 0.;
    double s = 0.;

    // single medium, spatially constant cross sections
    if (_config->hasSingleConstantMedium())
    {
        double section = state(0, 0).mix->sectionExt(pp->wavelength());
        while (generator->next())
        {
            if (generator->m() >= 0)
            {
                tau += section * state(generator->m(), 0).n * generator->ds();
                if (tau >= taumax) return std::numeric_limits<double>::infinity();
            }
            s += generator->ds();
            if (s > distance) break;
        }
    }

    // multiple media, spatially constant cross sections
    else if (_config->hasMultipleConstantMedia())
    {
        ShortArray<8> sectionv(_numMedia);
        for (int h = 0; h != _numMedia; ++h) sectionv[h] = state(0, h).mix->sectionExt(pp->wavelength());
        while (generator->next())
        {
            double ds = generator->ds();
            int m = generator->m();
            if (m >= 0)
            {
                for (int h = 0; h != _numMedia; ++h) tau += sectionv[h] * state(m, h).n * ds;
                if (tau >= taumax) return std::numeric_limits<double>::infinity();
            }
            s += ds;
            if (s > distance) break;
        }
    }

    // spatially variable cross sections
    else
    {
        while (generator->next())
        {
            double ds = generator->ds();
            int m = generator->m();
            if (m >= 0)
            {
                double lambda = pp->perceivedWavelength(state(m).v, _config->lyaExpansionRate() * s);
                tau += opacityExt(lambda, m) * ds;
                if (tau >= taumax) return std::numeric_limits<double>::infinity();
            }
            s += ds;
            if (s > distance) break;
        }
    }

    return tau;
}

////////////////////////////////////////////////////////////////////

void MediumSystem::clearRadiationField(bool primary)
{
    if (primary)
    {
        _rf1.setToZero();
        if (_rf2.size()) _rf2.setToZero();
    }
    else
    {
        _rf2c.setToZero();
    }
}

////////////////////////////////////////////////////////////////////

void MediumSystem::storeRadiationField(bool primary, int m, int ell, double Lds)
{
    if (primary)
        LockFree::add(_rf1(m, ell), Lds);
    else
        LockFree::add(_rf2c(m, ell), Lds);
}

////////////////////////////////////////////////////////////////////

void MediumSystem::communicateRadiationField(bool primary)
{
    if (primary)
        ProcessManager::sumToAll(_rf1.data());
    else
    {
        ProcessManager::sumToAll(_rf2c.data());
        _rf2 = _rf2c;
    }
}

////////////////////////////////////////////////////////////////////

double MediumSystem::radiationField(int m, int ell) const
{
    double rf = 0.;
    if (_rf1.size()) rf += _rf1(m, ell);
    if (_rf2.size()) rf += _rf2(m, ell);
    return rf;
}

////////////////////////////////////////////////////////////////////

Array MediumSystem::meanIntensity(int m) const
{
    int numWavelengths = _wavelengthGrid->numBins();
    Array Jv(numWavelengths);
    double factor = 1. / (4. * M_PI * volume(m));
    for (int ell = 0; ell < numWavelengths; ell++)
    {
        Jv[ell] = radiationField(m, ell) * factor / _wavelengthGrid->effectiveWidth(ell);
    }
    return Jv;
}

////////////////////////////////////////////////////////////////////

double MediumSystem::indicativeDustTemperature(int m) const
{
    const Array& Jv = meanIntensity(m);
    double sumRhoT = 0.;
    double sumRho = 0.;
    for (int h = 0; h != _numMedia; ++h)
    {
        if (isDust(h))
        {
            double rho = massDensity(m, h);
            if (rho > 0.)
            {
                double T = mix(m, h)->equilibriumTemperature(Jv);
                sumRhoT += rho * T;
                sumRho += rho;
            }
        }
    }
    if (sumRho > 0.)
        return sumRhoT / sumRho;
    else
        return 0.;
}

////////////////////////////////////////////////////////////////////

double MediumSystem::absorbedDustLuminosity(int m) const
{
    double Labs = 0.;
    int numWavelengths = _wavelengthGrid->numBins();
    for (int ell = 0; ell < numWavelengths; ell++)
    {
        Labs +=
            opacityAbs(_wavelengthGrid->wavelength(ell), m, MaterialMix::MaterialType::Dust) * radiationField(m, ell);
    }
    return Labs;
}

////////////////////////////////////////////////////////////////////

double MediumSystem::totalAbsorbedDustLuminosity(bool primary) const
{
    double Labs = 0.;
    int numWavelengths = _wavelengthGrid->numBins();
    for (int ell = 0; ell != numWavelengths; ++ell)
    {
        double lambda = _wavelengthGrid->wavelength(ell);
        for (int m = 0; m != _numCells; ++m)
        {
            double rf = primary ? _rf1(m, ell) : _rf2(m, ell);
            Labs += opacityAbs(lambda, m, MaterialMix::MaterialType::Dust) * rf;
        }
    }
    return Labs;
}

////////////////////////////////////////////////////////////////////
