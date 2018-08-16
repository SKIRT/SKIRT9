/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "MediumSystem.hpp"
#include "Configuration.hpp"
#include "FatalError.hpp"
#include "Log.hpp"
#include "MaterialMix.hpp"
#include "Parallel.hpp"
#include "ParallelFactory.hpp"
#include "ProcessManager.hpp"
#include "ShortArray.hpp"
#include "StringUtils.hpp"
#include "Table.hpp"

////////////////////////////////////////////////////////////////////

void MediumSystem::setupSelfAfter()
{
    SimulationItem::setupSelfAfter();
    auto log = find<Log>();
    auto parfac = find<ParallelFactory>();

    // accumulator for allocated memory in bytes
    size_t allocatedBytes = 0;


    // ----- volumes -----

    // allocate space for the grid cell volumes
    _numCells = _grid->numCells();
    if (_numCells<1) throw FATALERROR("The spatial grid must have at least one cell");
    _Vv.resize(_numCells);
    allocatedBytes += _numCells*sizeof(double);

    // calculate the grid cell volumes
    log->info("Calculating cell volumes...");
    parfac->parallelDistributed()->call(_numCells, [this](size_t firstIndex, size_t numIndices)
    {
        for (size_t m=firstIndex; m!=firstIndex+numIndices; ++m) _Vv[m] = _grid->volume(m);
    });
    ProcessManager::sumToAll(_Vv);

    // ----- cell states -----

    // if there are no media, do not create media states
    _numMedia = _media.size();
    if (_numMedia)
    {
        // allocate temporary space for the densities per cell and per medium
        // so that we can calculate them in parallel and communicate the result between processes
        Table<2> nvv(_numCells,_numMedia);      // indexed on m,h

        // calculate the density in the cells (by sampling in 100 random positions within the cell)
        log->info("Calculating densities for " + std::to_string(_numCells) + " cells...");
        log->infoSetElapsed();
        int numSamples = find<Configuration>()->numDensitySamples();
        parfac->parallelDistributed()->call(_numCells,
                                            [this, &nvv, log, numSamples](size_t firstIndex, size_t numIndices)
        {
            ShortArray<8> sumv(_numMedia);
            for (size_t m=firstIndex; m!=firstIndex+numIndices; ++m)
            {
                if (m%10000==0) log->infoIfElapsed("Calculated cell densities: ", m, _numCells);

                sumv.clear();
                for (int n=0; n<numSamples; n++)
                {
                    Position bfr = _grid->randomPositionInCell(m);
                    for (int h=0; h<_numMedia; h++) sumv[h] += _media[h]->numberDensity(bfr);
                }
                for (int h=0; h<_numMedia; h++) nvv(m,h) = sumv[h]/numSamples;
            }
        });
        ProcessManager::sumToAll(nvv.data());
        log->info("Done calculating cell densities");

        // allocate space for the state per cell and per medium
        _Svv.resize(_numCells*_numMedia);
        allocatedBytes += _numCells*_numMedia*sizeof(State);

        // insert the densities and the material mix pointers into the states
        for (int m=0; m!=_numCells; ++m)
        {
            Position bfr = _grid->centralPositionInCell(m);
            for (int h=0; h!=_numMedia; ++h)
            {
                state(m,h).n = nvv(m,h);
                state(m,h).mix = _media[h]->mix(bfr);
            }
        }
    }

    // ----------

    // log allocated memory size
    log->info(typeAndName() + " allocated " + StringUtils::toMemSizeString(allocatedBytes) + " of memory");
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
    return _media.size(); // don't use cached value because this function should work before setup has completed
}

////////////////////////////////////////////////////////////////////

int MediumSystem::numCells() const
{
    return _grid->numCells(); // don't use cached value because this function should work before setup has completed
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
