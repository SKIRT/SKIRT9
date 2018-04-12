/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "ParallelFactory.hpp"
#include "FatalError.hpp"

////////////////////////////////////////////////////////////////////

ParallelFactory::ParallelFactory()
{
}

////////////////////////////////////////////////////////////////////

ParallelFactory::ParallelFactory(SimulationItem* parent)
{
    parent->addChild(this);
}

////////////////////////////////////////////////////////////////////

void ParallelFactory::setMaxThreadCount(int value)
{
    _maxThreadCount = std::max(1, value);
}

////////////////////////////////////////////////////////////////////

int ParallelFactory::maxThreadCount() const
{
    return _maxThreadCount;
}

////////////////////////////////////////////////////////////////////

int ParallelFactory::defaultThreadCount()
{
    int count = std::thread::hardware_concurrency();
    return count>0 ? count : 1;
}

////////////////////////////////////////////////////////////////////

Parallel* ParallelFactory::parallel(int maxThreadCount)
{
    // Verify that we're being called from our parent thread
    if (std::this_thread::get_id() != _parentThread)
        throw FATALERROR("Parallel not spawned from thread that constructed the factory");

    // Determine the appropriate number of threads
    int numThreads = maxThreadCount>0 ? std::min(maxThreadCount, _maxThreadCount) : _maxThreadCount;

    // Get or create a child with that number of threads
    auto& child = _children[numThreads];
    if (!child) child.reset( new Parallel(numThreads, this) );
    return child.get();
}

////////////////////////////////////////////////////////////////////

int ParallelFactory::currentThreadIndex() const
{
    auto search = _indices.find(std::this_thread::get_id());
    if (search == _indices.end()) throw FATALERROR("Current thread index was not found");
    return search->second;
}

////////////////////////////////////////////////////////////////////

void ParallelFactory::addThreadIndex(std::thread::id threadId, int index)
{
    _indices[threadId] = index;
}

////////////////////////////////////////////////////////////////////
