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
#include "MaterialState.hpp"
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

    _numCells = _grid->numCells();
    if (_numCells < 1) throw FATALERROR("The spatial grid must have at least one cell");
    _numMedia = _media.size();
    size_t allocatedBytes = 0;

    // ----- allocate memory for the medium state -----

    // basic configuration
    _state.initConfiguration(_numCells, _numMedia);

    // common state variables
    vector<StateVariable> variables;
    variables.emplace_back(StateVariable::volume());
    if (_config->hasMovingMedia()) variables.emplace_back(StateVariable::bulkVelocity());
    if (_config->hasMagneticField()) variables.emplace_back(StateVariable::magneticField());
    _state.initCommonStateVariables(variables);

    // specific state variables
    for (auto medium : _media) _state.initSpecificStateVariables(medium->mix()->specificStateVariableInfo());

    // finalize
    allocatedBytes += _state.initAllocate() * sizeof(double);

    // ----- allocate memory for the radiation field -----

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

    // ----- cache info on the dust emission wavelength grid -----

    if (_config->hasDustEmission())
    {
        _numDustEmissionWavelengths = _config->dustEmissionWLG()->extlambdav().size();
    }

    // ----- obtain the material mix pointers -----

    if (_config->hasVariableMedia())
    {
        // we need a separate material mix pointer for each spatial cell (per medium component)
        _mixv.resize(_numCells * _numMedia);
        _mixPerCell = true;
        for (int m = 0; m != _numCells; ++m)
        {
            Position bfr = _grid->centralPositionInCell(m);
            for (int h = 0; h != _numMedia; ++h) _mixv[m * _numMedia + h] = _media[h]->mix(bfr);
        }
    }
    else
    {
        // the material mix pointer is identical for all spatial cells (per medium component)
        _mixv.resize(_numMedia);
        for (int h = 0; h != _numMedia; ++h) _mixv[h] = _media[h]->mix();
    }
    allocatedBytes += _mixv.size() * sizeof(MaterialMix*);

    // cache a list of medium component indices for each material type
    for (int h = 0; h != _numMedia; ++h)
    {
        switch (mix(0, h)->materialType())
        {
            case MaterialMix::MaterialType::Dust: _dust_hv.push_back(h); break;
            case MaterialMix::MaterialType::Gas: _gas_hv.push_back(h); break;
            case MaterialMix::MaterialType::Electrons: _elec_hv.push_back(h); break;
        }
    }

    // ----- inform user about allocated memory -----

    log->info(typeAndName() + " allocated " + StringUtils::toMemSizeString(allocatedBytes) + " of memory");

    // ----- calculate cell densities, bulk velocities, and volumes in parallel -----

    log->info("Calculating densities for " + std::to_string(_numCells) + " cells...");
    auto dic = _grid->interface<DensityInCellInterface>(0, false);  // optional fast-track interface for densities
    int numSamples = _config->numDensitySamples();
    log->infoSetElapsed(_numCells);
    parfac->parallelDistributed()->call(_numCells, [this, log, dic, numSamples](size_t firstIndex, size_t numIndices) {
        ShortArray nsumv(_numMedia);

        while (numIndices)
        {
            size_t currentChunkSize = min(logProgressChunkSize, numIndices);
            for (size_t m = firstIndex; m != firstIndex + currentChunkSize; ++m)
            {
                Position center = _grid->centralPositionInCell(m);

                // volume
                _state.setVolume(m, _grid->volume(m));

                // density: use optional fast-track interface or sample 100 random positions within the cell
                if (dic)
                {
                    for (int h = 0; h != _numMedia; ++h) _state.setNumberDensity(m, h, dic->numberDensity(h, m));
                }
                else
                {
                    nsumv.clear();
                    for (int n = 0; n < numSamples; n++)
                    {
                        Position bfr = _grid->randomPositionInCell(m);
                        for (int h = 0; h != _numMedia; ++h) nsumv[h] += _media[h]->numberDensity(bfr);
                    }
                    for (int h = 0; h != _numMedia; ++h) _state.setNumberDensity(m, h, nsumv[h] / numSamples);
                }

                // bulk velocity: weighted average at cell center; for oligochromatic simulations, leave at zero
                if (_config->hasMovingMedia())
                {
                    double n = 0.;
                    Vec v;
                    for (int h = 0; h != _numMedia; ++h)
                    {
                        n += _state.numberDensity(m, h);
                        v += _state.numberDensity(m, h) * _media[h]->bulkVelocity(center);
                    }
                    // leave bulk velocity at zero if cell has no material
                    if (n > 0.) _state.setBulkVelocity(m, v / n);
                }

                // magnetic field: retrieve from medium component that specifies it, if any
                if (_config->hasMagneticField())
                {
                    _state.setMagneticField(m, _media[_config->magneticFieldMediumIndex()]->magneticField(center));
                }

                // specific state variables other than density
                for (int h = 0; h != _numMedia; ++h)
                {
                    // retrieve input model temperature and parameters from medium component
                    MaterialState mst(_state, m, h);
                    double T = _media[h]->temperature(center);
                    Array params;
                    _media[h]->parameters(center, params);
                    mix(m, h)->initializeSpecificState(&mst, T, params);
                }
            }
            log->infoIfElapsed("Calculated cell densities: ", currentChunkSize);
            firstIndex += currentChunkSize;
            numIndices -= currentChunkSize;
        }
    });

    // communicate the calculated states across multiple processes, if needed
    _state.initCommunicate();

    log->info("Done calculating cell densities");
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

const MaterialMix* MediumSystem::mix(int m, int h) const
{
    return _mixPerCell ? _mixv[m * _numMedia + h] : _mixv[h];
}

////////////////////////////////////////////////////////////////////

bool MediumSystem::hasMaterialType(MaterialMix::MaterialType type) const
{
    for (int h = 0; h != _numMedia; ++h)
        if (mix(0, h)->materialType() == type) return true;
    return false;
}

////////////////////////////////////////////////////////////////////

bool MediumSystem::isMaterialType(MaterialMix::MaterialType type, int h) const
{
    return mix(0, h)->materialType() == type;
}

////////////////////////////////////////////////////////////////////

double MediumSystem::volume(int m) const
{
    return _state.volume(m);
}

////////////////////////////////////////////////////////////////////

Vec MediumSystem::bulkVelocity(int m) const
{
    return _state.bulkVelocity(m);
}

////////////////////////////////////////////////////////////////////

Vec MediumSystem::magneticField(int m) const
{
    return _state.magneticField(m);
}

////////////////////////////////////////////////////////////////////

double MediumSystem::numberDensity(int m, int h) const
{
    return _state.numberDensity(m, h);
}

////////////////////////////////////////////////////////////////////

double MediumSystem::massDensity(int m, int h) const
{
    return _state.numberDensity(m, h) * mix(m, h)->mass();
}

////////////////////////////////////////////////////////////////////

double MediumSystem::opacityAbs(double lambda, int m, int h) const
{
    MaterialState mst(_state, m, h);
    return mix(m, h)->opacityAbs(lambda, &mst, nullptr);
}

////////////////////////////////////////////////////////////////////

double MediumSystem::opacitySca(double lambda, int m, int h) const
{
    MaterialState mst(_state, m, h);
    return mix(m, h)->opacitySca(lambda, &mst, nullptr);
}

////////////////////////////////////////////////////////////////////

double MediumSystem::opacityExt(double lambda, int m, int h) const
{
    MaterialState mst(_state, m, h);
    return mix(m, h)->opacityExt(lambda, &mst, nullptr);
}

////////////////////////////////////////////////////////////////////

double MediumSystem::opacityAbs(double lambda, int m, MaterialMix::MaterialType type) const
{
    double result = 0.;
    for (int h = 0; h != _numMedia; ++h)
        if (mix(0, h)->materialType() == type) result += opacityAbs(lambda, m, h);
    return result;
}

////////////////////////////////////////////////////////////////////

double MediumSystem::opacityExt(double lambda, int m, MaterialMix::MaterialType type) const
{
    double result = 0.;
    for (int h = 0; h != _numMedia; ++h)
        if (mix(0, h)->materialType() == type) result += opacityExt(lambda, m, h);
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

double MediumSystem::perceivedWavelengthForScattering(const PhotonPacket* pp) const
{
    if (_config->hasMovingMedia())
        return pp->perceivedWavelength(_state.bulkVelocity(pp->interactionCellIndex()),
                                       _config->hubbleExpansionRate() * pp->interactionDistance());
    else
        return pp->wavelength();
}

////////////////////////////////////////////////////////////////////

double MediumSystem::albedoForScattering(const PhotonPacket* pp) const
{
    double lambda = perceivedWavelengthForScattering(pp);
    int m = pp->interactionCellIndex();
    if (m < 0) throw FATALERROR("Cannot locate photon packet interaction point");

    double ksca = 0.;
    double kext = 0.;
    for (int h = 0; h != _numMedia; ++h)
    {
        MaterialState mst(_state, m, h);
        ksca += mix(m, h)->opacitySca(lambda, &mst, pp);
        kext += mix(m, h)->opacityExt(lambda, &mst, pp);
    }
    return kext > 0. ? ksca / kext : 0.;
}

////////////////////////////////////////////////////////////////////

bool MediumSystem::weightsForScattering(ShortArray& wv, double lambda, const PhotonPacket* pp) const
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
        MaterialState mst(_state, m, h);
        wv[h] = mix(m, h)->opacitySca(lambda, &mst, pp);
        sum += wv[h];
    }

    // normalize the weights
    if (sum > 0.)
    {
        for (int h = 0; h != _numMedia; ++h) wv[h] /= sum;
        return true;
    }

    // none of the media scatter
    return false;
}

////////////////////////////////////////////////////////////////////

void MediumSystem::peelOffScattering(double lambda, const ShortArray& wv, Direction bfkobs, Direction bfky,
                                     PhotonPacket* pp, PhotonPacket* ppp) const
{
    // get the cell hosting the scattering event
    int m = pp->interactionCellIndex();

    // the outgoing wavelength (in the bulk velocity frame)
    double emissionLambda = lambda;

    // calculate the weighted sum of the effects on the Stokes vector and on the wavelength for all media
    double I = 0., Q = 0., U = 0., V = 0.;
    for (int h = 0; h != _numMedia; ++h)
    {
        double localLambda = lambda;
        MaterialState mst(_state, m, h);
        mix(m, h)->peeloffScattering(I, Q, U, V, localLambda, wv[h], bfkobs, bfky, &mst, pp);

        // if this material mix changed the wavelength, it is copied as the outgoing wavelength
        // if more than one material mix changes the wavelength, only the last one is preserved
        if (localLambda != lambda) emissionLambda = localLambda;
    }

    // pass the result to the peel-off photon packet
    ppp->launchScatteringPeelOff(pp, bfkobs, _state.bulkVelocity(m), emissionLambda, I);
    if (_config->hasPolarization()) ppp->setPolarized(I, Q, U, V, pp->normal());
}

////////////////////////////////////////////////////////////////////

void MediumSystem::simulateScattering(Random* random, PhotonPacket* pp) const
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
            MaterialState mst(_state, m, h);
            return mix(m, h)->opacitySca(lambda, &mst, pp);
        });

        // randomly select an index
        h = NR::locateClip(Xv, random->uniform());
    }

    // actually perform the scattering event for this cell and medium component
    MaterialState mst(_state, m, h);
    mix(m, h)->performScattering(lambda, &mst, pp);
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

double MediumSystem::getOpticalDepth(const SpatialGridPath* path, double lambda, MaterialMix::MaterialType type) const
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

void MediumSystem::setOpticalDepths(PhotonPacket* pp) const
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
    if (_config->hasSingleConstantSectionMedium())
    {
        double section = mix(0, 0)->sectionExt(pp->wavelength());
        for (auto& segment : pp->segments())
        {
            if (segment.m >= 0) tau += section * _state.numberDensity(segment.m, 0) * segment.ds;
            pp->setOpticalDepth(i++, tau);
        }
    }

    // multiple media, spatially constant cross sections
    else if (_config->hasMultipleConstantSectionMedia())
    {
        ShortArray sectionv(_numMedia);
        for (int h = 0; h != _numMedia; ++h) sectionv[h] = mix(0, h)->sectionExt(pp->wavelength());
        for (auto& segment : pp->segments())
        {
            if (segment.m >= 0)
                for (int h = 0; h != _numMedia; ++h)
                    tau += sectionv[h] * _state.numberDensity(segment.m, h) * segment.ds;
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
                double lambda =
                    pp->perceivedWavelength(_state.bulkVelocity(segment.m), _config->hubbleExpansionRate() * segment.s);
                tau += opacityExt(lambda, segment.m) * segment.ds;
            }
            pp->setOpticalDepth(i++, tau);
        }
    }
}

////////////////////////////////////////////////////////////////////

bool MediumSystem::setInteractionPoint(PhotonPacket* pp, double tauscat) const
{
    auto generator = getPathSegmentGenerator(_grid, pp);
    double tau = 0.;
    double s = 0.;

    // loop over the segments of the path until the interaction optical depth is reached or the path ends

    // --> single medium, spatially constant cross sections
    if (_config->hasSingleConstantSectionMedium())
    {
        double section = mix(0, 0)->sectionExt(pp->wavelength());
        while (generator->next())
        {
            // remember the cumulative optical depth and distance at the start of this segment
            // so that we can interpolate the interaction point should it happen to be inside this segment
            double tau0 = tau;
            double s0 = s;

            // calculate the cumulative optical depth and distance at the end of this segment
            double ds = generator->ds();
            int m = generator->m();
            if (m >= 0) tau += section * _state.numberDensity(m, 0) * generator->ds();
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
    else if (_config->hasMultipleConstantSectionMedia())
    {
        ShortArray sectionv(_numMedia);
        for (int h = 0; h != _numMedia; ++h) sectionv[h] = mix(0, h)->sectionExt(pp->wavelength());
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
                for (int h = 0; h != _numMedia; ++h) tau += sectionv[h] * _state.numberDensity(m, h) * ds;
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
                double lambda = pp->perceivedWavelength(_state.bulkVelocity(m), _config->hubbleExpansionRate() * s);
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

double MediumSystem::getOpticalDepth(PhotonPacket* pp, double distance) const
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
    if (_config->hasSingleConstantSectionMedium())
    {
        double section = mix(0, 0)->sectionExt(pp->wavelength());
        while (generator->next())
        {
            if (generator->m() >= 0)
            {
                tau += section * _state.numberDensity(generator->m(), 0) * generator->ds();
                if (tau >= taumax) return std::numeric_limits<double>::infinity();
            }
            s += generator->ds();
            if (s > distance) break;
        }
    }

    // multiple media, spatially constant cross sections
    else if (_config->hasMultipleConstantSectionMedia())
    {
        ShortArray sectionv(_numMedia);
        for (int h = 0; h != _numMedia; ++h) sectionv[h] = mix(0, h)->sectionExt(pp->wavelength());
        while (generator->next())
        {
            double ds = generator->ds();
            int m = generator->m();
            if (m >= 0)
            {
                for (int h = 0; h != _numMedia; ++h) tau += sectionv[h] * _state.numberDensity(m, h) * ds;
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
                double lambda = pp->perceivedWavelength(_state.bulkVelocity(m), _config->hubbleExpansionRate() * s);
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
    double factor = 1. / (4. * M_PI * _state.volume(m));
    for (int ell = 0; ell < numWavelengths; ell++)
    {
        Jv[ell] = radiationField(m, ell) * factor / _wavelengthGrid->effectiveWidth(ell);
    }
    return Jv;
}

////////////////////////////////////////////////////////////////////

double MediumSystem::indicativeTemperature(int m, MaterialMix::MaterialType type) const
{
    // get the radiation field, if available
    Array Jv;
    if (_config->hasRadiationField()) Jv = meanIntensity(m);

    // obtain the temperature and weight for each component of the requested type
    double sumRhoT = 0.;
    double sumRho = 0.;
    for (int h = 0; h != _numMedia; ++h)
    {
        if (isMaterialType(type, h))
        {
            double rho = massDensity(m, h);
            if (rho > 0.)
            {
                MaterialState mst(_state, m, h);
                double T = mix(m, h)->indicativeTemperature(&mst, Jv);
                sumRhoT += rho * T;
                sumRho += rho;
            }
        }
    }

    // return the average, if there is one
    return sumRho > 0. ? sumRhoT / sumRho : 0.;
}

////////////////////////////////////////////////////////////////////

double MediumSystem::indicativeDustTemperature(int m) const
{
    return indicativeTemperature(m, MaterialMix::MaterialType::Dust);
}

////////////////////////////////////////////////////////////////////

double MediumSystem::indicativeGasTemperature(int m) const
{
    return indicativeTemperature(m, MaterialMix::MaterialType::Gas);
}

////////////////////////////////////////////////////////////////////

double MediumSystem::dustLuminosity(int m) const
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

Array MediumSystem::dustEmissionSpectrum(int m) const
{
    const Array& Jv = meanIntensity(m);
    Array ev(_numDustEmissionWavelengths);
    for (int h : _dust_hv)
    {
        MaterialState mst(_state, m, h);
        ev += mix(m, h)->emissionSpectrum(&mst, Jv);
    }
    return ev;
}

////////////////////////////////////////////////////////////////////
