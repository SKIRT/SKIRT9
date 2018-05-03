/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "MultiProcessParallel.hpp"
#include "FatalError.hpp"
#include "ProcessManager.hpp"

////////////////////////////////////////////////////////////////////

MultiProcessParallel::MultiProcessParallel(int /*threadCount*/)
{
    if (ProcessManager::isRoot())
    {
        constructThreads(1);
    }
}

////////////////////////////////////////////////////////////////////

MultiProcessParallel::~MultiProcessParallel()
{
    if (ProcessManager::isRoot())
    {
        destroyThreads();
    }
}

////////////////////////////////////////////////////////////////////

void MultiProcessParallel::call(std::function<void(size_t,size_t)> target, size_t maxIndex)
{
    // In the root process, the child thread performs work, and the parent thread serves other processes
    if (ProcessManager::isRoot())
    {
        // Copy the target function so it can be invoked from the child thread
        _target = target;

        // Initialize the chunk maker
        _chunkMaker.initialize(maxIndex, 1, ProcessManager::size());

        // Activate child thread
        activateThreads();

        // Serve chunks to other processes
        int rank = 0;
        while (true)
        {
            rank = ProcessManager::waitForChunkRequest();
            size_t firstIndex, numIndices;
            if (!_chunkMaker.next(firstIndex, numIndices)) break;
            ProcessManager::serveChunkRequest(rank, firstIndex, numIndices);
        }

        // Serve each non-root process an empty chunk as a terminating signal
        ProcessManager::serveChunkRequest(rank, 0, 0);
        for (int i = 2; i!=ProcessManager::size(); ++i)
        {
            rank = ProcessManager::waitForChunkRequest();
            ProcessManager::serveChunkRequest(rank, 0, 0);
        }

        // wait for our child thread to finish as well
        waitForThreads();
    }

    // In non-root processes, the parent (and only) thread performs work in a straightforward loop
    else
    {
        size_t firstIndex, numIndices;
        while (ProcessManager::requestChunk(firstIndex, numIndices)) target(firstIndex, numIndices);
    }
}

////////////////////////////////////////////////////////////////////

bool MultiProcessParallel::doSomeWork()
{
    return _chunkMaker.callForNext(_target);
}

///////////////////////////////////////////////////////////////////
