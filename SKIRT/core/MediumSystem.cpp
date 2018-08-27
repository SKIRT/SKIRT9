/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "MediumSystem.hpp"
#include "Configuration.hpp"
#include "FatalError.hpp"
#include "Log.hpp"
#include "MaterialMix.hpp"
#include "NR.hpp"
#include "Parallel.hpp"
#include "ParallelFactory.hpp"
#include "PhotonPacket.hpp"
#include "ProcessManager.hpp"
#include "Random.hpp"
#include "ShortArray.hpp"
#include "StringUtils.hpp"
#include "Table.hpp"

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

    // ----- allocate memory -----

    _numCells = _grid->numCells();
    if (_numCells<1) throw FATALERROR("The spatial grid must have at least one cell");
    _numMedia = _media.size();

    size_t allocatedBytes = 0;
    _state1v.resize(_numCells);
    allocatedBytes += _state1v.size()*sizeof(State1);
    _state2vv.resize(_numCells*_numMedia);
    allocatedBytes += _state2vv.size()*sizeof(State2);
    log->info(typeAndName() + " allocated " + StringUtils::toMemSizeString(allocatedBytes) + " of memory");

    // ----- calculate cell densities, bulk velocities, and volumes in parallel -----

    log->info("Calculating densities for " + std::to_string(_numCells) + " cells...");
    int numSamples = find<Configuration>()->numDensitySamples();
    log->infoSetElapsed(_numCells);
    parfac->parallelDistributed()->call(_numCells,
                                        [this, log, numSamples](size_t firstIndex, size_t numIndices)
    {
        ShortArray<8> nsumv(_numMedia);

        while (numIndices)
        {
            size_t currentChunkSize = min(logProgressChunkSize, numIndices);
            for (size_t m=firstIndex; m!=firstIndex+currentChunkSize; ++m)
            {
                // density: sample 100 random positions within the cell
                nsumv.clear();
                for (int n=0; n<numSamples; n++)
                {
                    Position bfr = _grid->randomPositionInCell(m);
                    for (int h=0; h<_numMedia; h++) nsumv[h] += _media[h]->numberDensity(bfr);
                }
                for (int h=0; h<_numMedia; h++) state(m,h).n = nsumv[h]/numSamples;

                // bulk velocity: weighted average; assumes densities have been calculated
                Position bfr = _grid->centralPositionInCell(m);
                double n = 0.;
                Vec v;
                for (int h=0; h!=_numMedia; ++h)
                {
                    n += state(m,h).n;
                    v += state(m,h).n * _media[h]->bulkVelocity(bfr);
                }
                if (n > 0.) state(m).v = v / n;  // leave bulk velocity at zero if cell has no material

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

    for (int m=0; m!=_numCells; ++m)
    {
        Position bfr = _grid->centralPositionInCell(m);
        for (int h=0; h!=_numMedia; ++h) state(m,h).mix = _media[h]->mix(bfr);
    }
}

////////////////////////////////////////////////////////////////////

void MediumSystem::communicateStates()
{
    if (!ProcessManager::isMultiProc()) return;

    // NOTE: once the design of the state data structures is stable, a custom communication procedure could be provided
    //       in the meantime, we copy the data into a temporary table so we can use the standard sumToAll procedure
    Table<2> data;

    // volumes and bulk velocities
    data.resize(_numCells,4);
    for (int m=0; m!=_numCells; ++m)
    {
        data(m,0) = state(m).V;
        data(m,1) = state(m).v.x();
        data(m,2) = state(m).v.y();
        data(m,3) = state(m).v.z();
    }
    ProcessManager::sumToAll(data.data());
    for (int m=0; m!=_numCells; ++m)
    {
        state(m).V = data(m,0);
        state(m).v = Vec(data(m,1), data(m,3), data(m,3));
    }

    // densities
    data.resize(_numCells,_numMedia);
    for (int m=0; m!=_numCells; ++m) for (int h=0; h!=_numMedia; ++h) data(m,h) = state(m,h).n;
    ProcessManager::sumToAll(data.data());
    for (int m=0; m!=_numCells; ++m) for (int h=0; h!=_numMedia; ++h) state(m,h).n = data(m,h);
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

bool MediumSystem::hasMaterialType(MaterialMix::MaterialType type) const
{
    for (int h=0; h!=_numMedia; ++h) if (state(0,h).mix->materialType() == type) return true;
    return false;
}

////////////////////////////////////////////////////////////////////

bool MediumSystem::isMaterialType(MaterialMix::MaterialType type, int h) const
{
    return state(0,h).mix->materialType() == type;
}

////////////////////////////////////////////////////////////////////

double MediumSystem::numberDensity(int m, int h) const
{
    return state(m,h).n;
}

////////////////////////////////////////////////////////////////////

double MediumSystem::massDensity(int m, int h) const
{
    return state(m,h).n * state(m,h).mix->mass();
}

////////////////////////////////////////////////////////////////////

const MaterialMix* MediumSystem::mix(int m, int h) const
{
    return state(m,h).mix;
}

////////////////////////////////////////////////////////////////////

const MaterialMix* MediumSystem::randomMixForScattering(Random* random, double lambda, int m) const
{
    int h = 0;
    if (_numMedia>1)
    {
        Array Xv;
        NR::cdf(Xv, _numMedia, [this,lambda,m](int h){ return state(m,h).n * state(m,h).mix->sectionSca(lambda); });
        h = NR::locateClip(Xv, random->uniform());
    }
    return state(m,h).mix;
}

////////////////////////////////////////////////////////////////////

double MediumSystem::opacitySca(double lambda, int m, int h) const
{
    return state(m,h).n * state(m,h).mix->sectionSca(lambda);
}

////////////////////////////////////////////////////////////////////

double MediumSystem::opacitySca(double lambda, int m) const
{
    double result = 0.;
    for (int h=0; h!=_numMedia; ++h) result += state(m,h).n * state(m,h).mix->sectionSca(lambda);
    return result;
}

////////////////////////////////////////////////////////////////////

double MediumSystem::opacityExt(double lambda, int m, int h) const
{
    return state(m,h).n * state(m,h).mix->sectionExt(lambda);
}

////////////////////////////////////////////////////////////////////

double MediumSystem::opacityExt(double lambda, int m) const
{
    double result = 0.;
    for (int h=0; h!=_numMedia; ++h) result += state(m,h).n * state(m,h).mix->sectionExt(lambda);
    return result;
}

////////////////////////////////////////////////////////////////////

double MediumSystem::opacityExt(double lambda, int m, MaterialMix::MaterialType type) const
{
    double result = 0.;
    for (int h=0; h!=_numMedia; ++h)
        if (state(0,h).mix->materialType() == type) result += state(m,h).n * state(m,h).mix->sectionExt(lambda);
    return result;
}

////////////////////////////////////////////////////////////////////

double MediumSystem::albedo(double lambda, int m, int h) const
{
    return state(m,h).mix->albedo(lambda);
}

////////////////////////////////////////////////////////////////////

double MediumSystem::albedo(double lambda, int m) const
{
    double ksca = 0.;
    double kext = 0.;
    for (int h=0; h!=_numMedia; ++h)
    {
        double n = state(m,h).n;
        auto mix = state(m,h).mix;
        ksca += n * mix->sectionSca(lambda);
        kext += n * mix->sectionExt(lambda);
    }
    return kext>0. ? ksca/kext : 0.;
}

////////////////////////////////////////////////////////////////////

double MediumSystem::opticalDepth(PhotonPacket* pp, double distance)
{
    // determine the path and store the geometric details in the photon packet
    _grid->path(pp);

    // calculate and return the optical depth at the specified distance
    return pp->opticalDepth([this,pp](int m){ return opacityExt(pp->perceivedWavelength(state(m).v), m); }, distance);
}

////////////////////////////////////////////////////////////////////

void MediumSystem::fillOpticalDepth(PhotonPacket* pp)
{
    // determine the path and store the geometric details in the photon packet
    _grid->path(pp);

    // calculate and store the optical depth details in the photon package
    pp->fillOpticalDepth([this,pp](int m){ return opacityExt(pp->perceivedWavelength(state(m).v), m); });

    // verify that the result makes sense
    double tau = pp->tau();
    if (tau<0. || !std::isfinite(tau))
        throw FATALERROR("The optical depth along the path is not a positive number: tau = "
                         + StringUtils::toString(tau));
}

////////////////////////////////////////////////////////////////////
