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
#include "NR.hpp"
#include "Parallel.hpp"
#include "ParallelFactory.hpp"
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

Vec MediumSystem::bulkVelocity(int m)
{
    return state(m).v;
}

////////////////////////////////////////////////////////////////////

Vec MediumSystem::magneticField(int m)
{
    return state(m).B;
}

////////////////////////////////////////////////////////////////////

double MediumSystem::gasTemperature(int m)
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

const MaterialMix* MediumSystem::randomMixForScattering(Random* random, double lambda, int m) const
{
    int h = 0;
    if (_numMedia > 1)
    {
        Array Xv;
        NR::cdf(Xv, _numMedia, [this, lambda, m](int h) { return opacitySca(lambda, m, h); });
        h = NR::locateClip(Xv, random->uniform());
    }
    return state(m, h).mix;
}

////////////////////////////////////////////////////////////////////

double MediumSystem::opacitySca(double lambda, int m, int h) const
{
    double n = state(m, h).n;
    if (n <= 0.)
        return 0.;
    else if (h == _config->lyaMediumIndex())
        return n * LyaUtils::section(lambda, state(m).T);
    else
        return n * state(m, h).mix->sectionSca(lambda);
}

////////////////////////////////////////////////////////////////////

double MediumSystem::opacityAbs(double lambda, int m, int h) const
{
    // no need to check for Lyman-alpha because the material mix returns the correct zero value
    double n = state(m, h).n;
    if (n <= 0.)
        return 0.;
    else
        return n * state(m, h).mix->sectionAbs(lambda);
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

double MediumSystem::opacityExt(double lambda, int m, int h) const
{
    double n = state(m, h).n;
    if (n <= 0.)
        return 0.;
    else if (h == _config->lyaMediumIndex())
        return n * LyaUtils::section(lambda, state(m).T);
    else
        return n * state(m, h).mix->sectionExt(lambda);
}

////////////////////////////////////////////////////////////////////

double MediumSystem::opacityExt(double lambda, int m) const
{
    double result = 0.;
    for (int h = 0; h != _numMedia; ++h) result += opacityExt(lambda, m, h);
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

double MediumSystem::opticalDepth(SpatialGridPath* path, double lambda, MaterialMix::MaterialType type)
{
    // determine the geometric details of the path
    _grid->path(path);

    // calculate the optical depth
    double tau = 0.;
    for (const auto& segment : path->segments())
    {
        if (segment.m >= 0) tau += opacityExt(lambda, segment.m, type) * segment.ds;
    }
    return tau;
}

////////////////////////////////////////////////////////////////////

namespace
{
    // The dynamic range of a positive IEEE 754 double precision floating point number is just smaller than 2^2098.
    // Consequently, the intensity of a photon packet will always become numerically zero beyond the point where the
    // cumulative optical depth along the path reaches tau = 2098 ln 2.
    constexpr double TAU_MAX = 2098 * M_LN2;
}

////////////////////////////////////////////////////////////////////

void MediumSystem::opticalDepth(PhotonPacket* pp)
{
    // determine the geometric details of the path
    _grid->path(pp);

    // calculate the cumulative optical depth and store it in the photon packet for each path segment;
    // because this function is at the heart of the photon life cycle, we implement various optimized versions
    double tau = 0.;
    int i = 0;

    // no kinematics and material properties are spatially constant
    if (!_config->hasMovingMedia() && !_config->hasVariableMedia())
    {
        // single medium (no kinematics, spatially constant)
        if (_numMedia == 1)
        {
            double section = state(0, 0).mix->sectionExt(pp->wavelength());
            for (auto& segment : pp->segments())
            {
                if (segment.m >= 0) tau += section * state(segment.m, 0).n * segment.ds;
                pp->setOpticalDepth(i++, tau);
            }
        }
        // multiple media (no kinematics, spatially constant)
        else
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
    }
    // with kinematics and/or spatially variable material properties
    else
    {
        for (auto& segment : pp->segments())
        {
            if (segment.m >= 0)
            {
                double lambda = pp->perceivedWavelength(state(segment.m).v, _config->lyaExpansionRate() * segment.s);
                tau += opacityExt(lambda, segment.m) * segment.ds;
                if (tau >= TAU_MAX)
                {
                    pp->setTerminalOpticalDepth(i, tau);
                    break;
                }
            }
            pp->setOpticalDepth(i++, tau);
        }
    }
}

////////////////////////////////////////////////////////////////////

double MediumSystem::opticalDepth(PhotonPacket* pp, double distance)
{
    // determine the geometric details of the path
    _grid->path(pp);

    // calculate the cumulative optical depth
    // because this function is at the heart of the photon life cycle, we implement various optimized versions
    double tau = 0.;

    // no kinematics and material properties are spatially constant
    if (!_config->hasMovingMedia() && !_config->hasVariableMedia())
    {
        // single medium (no kinematics, spatially constant)
        if (_numMedia == 1)
        {
            double section = state(0, 0).mix->sectionExt(pp->wavelength());
            for (auto& segment : pp->segments())
            {
                if (segment.m >= 0) tau += section * state(segment.m, 0).n * segment.ds;
                if (segment.s > distance) break;
            }
        }
        // multiple media (no kinematics, spatially constant)
        else
        {
            ShortArray<8> sectionv(_numMedia);
            for (int h = 0; h != _numMedia; ++h) sectionv[h] = state(0, h).mix->sectionExt(pp->wavelength());
            for (auto& segment : pp->segments())
            {
                if (segment.m >= 0)
                    for (int h = 0; h != _numMedia; ++h) tau += sectionv[h] * state(segment.m, h).n * segment.ds;
                if (segment.s > distance) break;
            }
        }
    }
    // with kinematics and/or spatially variable material properties
    else
    {
        for (auto& segment : pp->segments())
        {
            if (segment.m >= 0)
            {
                double lambda = pp->perceivedWavelength(state(segment.m).v, _config->lyaExpansionRate() * segment.s);
                tau += opacityExt(lambda, segment.m) * segment.ds;
                if (tau >= TAU_MAX) break;
            }
            if (segment.s > distance) break;
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

double MediumSystem::totalAbsorbedLuminosity(bool primary, MaterialMix::MaterialType type) const
{
    double Labs = 0.;
    int numWavelengths = _wavelengthGrid->numBins();
    for (int ell = 0; ell != numWavelengths; ++ell)
    {
        double lambda = _wavelengthGrid->wavelength(ell);
        for (int m = 0; m != _numCells; ++m)
        {
            double rf = primary ? _rf1(m, ell) : _rf2(m, ell);
            Labs += opacityAbs(lambda, m, type) * rf;
        }
    }
    return Labs;
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

double MediumSystem::absorbedLuminosity(int m, MaterialMix::MaterialType type) const
{
    double Labs = 0.;
    int numWavelengths = _wavelengthGrid->numBins();
    for (int ell = 0; ell < numWavelengths; ell++)
    {
        Labs += opacityAbs(_wavelengthGrid->wavelength(ell), m, type) * radiationField(m, ell);
    }
    return Labs;
}

////////////////////////////////////////////////////////////////////
