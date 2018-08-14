/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "MediumSystem.hpp"
#include "Configuration.hpp"
#include "Log.hpp"
#include "Parallel.hpp"
#include "ParallelFactory.hpp"
#include "ProcessManager.hpp"
#include "ShortArray.hpp"
#include "StringUtils.hpp"

////////////////////////////////////////////////////////////////////

void MediumSystem::setupSelfAfter()
{
    SimulationItem::setupSelfAfter();
    auto log = find<Log>();
    auto parfac = find<ParallelFactory>();

    // allocate the tables that hold essential dust cell properties
    int numCells = _grid->numCells();
    int numMedia = _media.size();
    _Vv.resize(numCells);
    _nvv.resize(numCells,numMedia);

    // calculate and log allocated memory size
    size_t memSize = 0;
    memSize += _Vv.size();
    memSize += _nvv.size();
    log->info(typeAndName() + " allocated " + StringUtils::toMemSizeString(memSize*sizeof(double)) + " of memory");

    // calculate the volume of the cells
    log->info("Calculating cell volumes...");
    parfac->parallelDistributed()->call(numCells, [this](size_t firstIndex, size_t numIndices)
    {
        for (size_t m=firstIndex; m!=firstIndex+numIndices; ++m) _Vv[m] = _grid->volume(m);
    });
    ProcessManager::sumToAll(_Vv);

    // calculate the density in the cells (by sampling in 100 random positions within the cell)
    log->info("Calculating densities for " + std::to_string(numCells) + " cells...");
    log->infoSetElapsed();
    int numSamples = find<Configuration>()->numDensitySamples();
    parfac->parallelDistributed()->call(numCells,
        [this, log, numCells, numMedia, numSamples](size_t firstIndex, size_t numIndices)
    {
        ShortArray<8> sumv(numMedia);
        for (size_t m=firstIndex; m!=firstIndex+numIndices; ++m)
        {
            if (m%10000==0) log->infoIfElapsed("Calculated cell densities: ", m, numCells);

            sumv.clear();
            for (int n=0; n<numSamples; n++)
            {
                Position bfr = _grid->randomPositionInCell(m);
                for (int h=0; h<numMedia; h++) sumv[h] += _media[h]->numberDensity(bfr);
            }
            for (int h=0; h<numMedia; h++) _nvv(m,h) = sumv[h]/numSamples;
        }
    });
    ProcessManager::sumToAll(_nvv.data());
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
    return _media.size();
}

////////////////////////////////////////////////////////////////////

int MediumSystem::numCells() const
{
    return _Vv.size();
}

////////////////////////////////////////////////////////////////////
