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

void MultiThreadParallel::call(size_t maxIndex, std::function<void(size_t, size_t)> target)
{
    // Copy the target function so it can be invoked from any of the threads
    _target = target;

    // Initialize the chunk maker
    _chunkMaker.initialize(maxIndex, numThreads());

    // Activate child threads and wait until they are done; we don't do anything in the parent thread
    activateThreads();
    waitForThreads();
}

////////////////////////////////////////////////////////////////////

bool MultiThreadParallel::doSomeWork()
{
    return _chunkMaker.callForNext(_target);
}

////////////////////////////////////////////////////////////////////
