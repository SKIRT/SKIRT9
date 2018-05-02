/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "MultiThreadParallel.hpp"

////////////////////////////////////////////////////////////////////

MultiThreadParallel::MultiThreadParallel(int threadCount)
{
    constructThreads(threadCount);
}

////////////////////////////////////////////////////////////////////

MultiThreadParallel::~MultiThreadParallel()
{
    destroyThreads();
}

////////////////////////////////////////////////////////////////////

void MultiThreadParallel::call(std::function<void(size_t,size_t)> target, size_t maxIndex)
{
    // Copy the target function so it can be invoked from any of the threads
    _target = target;

    // Determine the chunk size
    const size_t numChunksPerThread = 8;   // empirical multiplicator to achieve acceptable load balancing
    _chunkSize = max(static_cast<size_t>(1), maxIndex / (numThreads()*numChunksPerThread));

    // Initialize the other data members
    _maxIndex = maxIndex;
    _nextIndex = 0;

    // Activate child threads and wait until they are done; we don't do anything in the parent thread
    activateThreads();
    waitForThreads();
}

////////////////////////////////////////////////////////////////////

bool MultiThreadParallel::doSomeWork()
{
    // Get and increment the next chunk index atomically, and report "done" if no more chunks are available
    size_t firstIndex = _nextIndex.fetch_add(_chunkSize);
    if (firstIndex >= _maxIndex) return false;

    // Invoke the target function
    _target(firstIndex, min(_chunkSize,_maxIndex-firstIndex));
    return true;
}

////////////////////////////////////////////////////////////////////
