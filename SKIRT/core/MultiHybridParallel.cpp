/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "MultiHybridParallel.hpp"
#include "FatalError.hpp"
#include "ProcessManager.hpp"

////////////////////////////////////////////////////////////////////

MultiHybridParallel::MultiHybridParallel(int threadCount)
{
    constructThreads(threadCount);
}

////////////////////////////////////////////////////////////////////

MultiHybridParallel::~MultiHybridParallel()
{
    destroyThreads();
}

////////////////////////////////////////////////////////////////////

void MultiHybridParallel::call(size_t maxIndex, std::function<void(size_t, size_t)> target)
{
    // Copy the target function so it can be invoked from the child threads
    _target = target;

    // In the root process, the parent thread serves chunks to other processes
    if (ProcessManager::isRoot())
    {
        // Initialize the chunk maker
        _chunkMaker.initialize(maxIndex, numThreads(), ProcessManager::size());

        // Activate child threads
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
        for (int i = 2; i != ProcessManager::size(); ++i)
        {
            rank = ProcessManager::waitForChunkRequest();
            ProcessManager::serveChunkRequest(rank, 0, 0);
        }

        // wait for our child threads to finish as well
        waitForThreads();
    }

    // In non-root processes, the parent thread requests chunks from the root process
    else
    {
        // Initialize the variables used to synchronize chunk requests with the child threads
        _requests = 0;
        _ready = false;
        _done = false;

        // Activate child threads
        activateThreads();

        // Serve chunks to the child threads upon their request
        bool success = true;
        while (success)
        {
            // Wait for new chunk request
            {
                std::unique_lock<std::mutex> lock(_mutex);
                while (!_requests || _ready) _conditionParent.wait(lock);
            }

            // Request a new chunk from the root process
            size_t firstIndex, numIndices;
            success = ProcessManager::requestChunk(firstIndex, numIndices);

            // Serve the chunk to one of our child threads, or tell our child threads that there are no more chunks
            {
                std::unique_lock<std::mutex> lock(_mutex);
                _firstIndex = firstIndex;
                _numIndices = numIndices;
                if (success)
                    _ready = true;
                else
                    _done = true;
            }
            _conditionChildren.notify_all();
        }

        // wait for our child threads to finish as well
        waitForThreads();
    }
}

////////////////////////////////////////////////////////////////////

bool MultiHybridParallel::doSomeWork()
{
    // In the root process, we share the chunk maker with the parent thread
    if (ProcessManager::isRoot())
    {
        return _chunkMaker.callForNext(_target);
    }

    // In non-root processes, we ask the parent thread to request a chunk from the root process
    else
    {
        // Request a new chunk
        {
            std::unique_lock<std::mutex> lock(_mutex);
            if (_done) return false;  // exit if there are no more chunks
            _requests++;
        }
        _conditionParent.notify_all();

        // Get the new chunk
        size_t firstIndex, numIndices;
        {
            std::unique_lock<std::mutex> lock(_mutex);
            while (!_ready && !_done) _conditionChildren.wait(lock);

            if (_done) return false;  // exit if there are no more chunks
            firstIndex = _firstIndex;
            numIndices = _numIndices;
            _ready = false;
            _requests--;
        }
        _conditionParent.notify_all();

        // Invoke the target function
        _target(firstIndex, numIndices);
        return true;
    }
}

///////////////////////////////////////////////////////////////////
