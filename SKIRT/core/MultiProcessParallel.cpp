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
        // Launch the child thread in a critical section
        {
            std::unique_lock<std::mutex> lock(_mutex);
            _childThread = std::thread(&MultiProcessParallel::run, this);
        }

        // Wait until the child thread is ready
        waitForChildThread();
    }
}

////////////////////////////////////////////////////////////////////

MultiProcessParallel::~MultiProcessParallel()
{
    if (ProcessManager::isRoot())
    {
        // Ask the parallel thread to exit in a critical section
        {
            std::unique_lock<std::mutex> lock(_mutex);
            _terminate = true;
            _conditionChild.notify_all();
        }

        // Wait for it to do so
         _childThread.join();
    }
}

////////////////////////////////////////////////////////////////////

void MultiProcessParallel::call(std::function<void(size_t,size_t)> target, size_t maxIndex)
{
    // In the root process, the child thread performs work, and the parent thread serves other processes
    if (ProcessManager::isRoot())
    {
        // Initialize shared data members and activate threads in a critical section
        {
            std::unique_lock<std::mutex> lock(_mutex);

            // Copy the target function so it can be invoked from the child thread
            _target = target;

            // Determine the chunk size
            const size_t numChunksPerThread = 8;   // empirical multiplicator to achieve acceptable load balancing
            _chunkSize = max(static_cast<size_t>(1), maxIndex / (ProcessManager::size()*numChunksPerThread));

            // Initialize the other data members
            _maxIndex = maxIndex;
            _exception = nullptr;
            _active = true;
            _nextIndex = 0;

            // Wake the child thread
            _conditionChild.notify_all();
        }

        // Serve chunks to other processes
        int rank = 0;
        while (true)
        {
            rank = ProcessManager::waitForChunkRequest();
            size_t firstIndex = _nextIndex.fetch_add(_chunkSize); // get and increment the next chunk index atomically
            if (firstIndex >= _maxIndex) break;                   // break if no more chunks are available
            ProcessManager::serveChunkRequest(rank, firstIndex, min(_chunkSize,_maxIndex-firstIndex));
        }

        // Serve each non-root process an empty chunk as a terminating signal
        ProcessManager::serveChunkRequest(rank, 0, 0);
        for (int i = 2; i!=ProcessManager::size(); ++i)
        {
            rank = ProcessManager::waitForChunkRequest();
            ProcessManager::serveChunkRequest(rank, 0, 0);
        }

        // Check for and process any exception raised in the child process
        if (_exception)
        {
            throw *_exception;  // throw by value (the memory for the heap-allocated exception is leaked)
        }
    }

    // In non-root processes, the parent (and only) thread performs work in a straightforward loop
    else
    {
        size_t firstIndex, numIndices;
        while (ProcessManager::requestChunk(firstIndex, numIndices)) target(firstIndex, numIndices);
    }
}

////////////////////////////////////////////////////////////////////

void MultiProcessParallel::run()
{
    while (true)
    {
        // Wait for new work in a critical section
        {
            std::unique_lock<std::mutex> lock(_mutex);

            // Indicate that this thread is no longer doing work
            _active = false;

            // Tell the main thread that the child thread is inactive
            _conditionMain.notify_all();

            // Wait for new work
            while (true)
            {
                _conditionChild.wait(lock);

                // Check for termination request (don't bother with _active; it's no longer used)
                if (_terminate) return;

                // Check that we actually have new work
                if (_active) break;
            }
        }

        // Do work as long as some is available
        try
        {
            while (true)
            {
                // Get and increment the next chunk index atomically, and break if no more chunks are available
                size_t firstIndex = _nextIndex.fetch_add(_chunkSize);
                if (firstIndex >= _maxIndex) break;

                // Invoke the target function
                _target(firstIndex, min(_chunkSize,_maxIndex-firstIndex));
            }
        }
        catch (FatalError& error)
        {
            // Make a copy of the exception
            reportException(new FatalError(error));
        }
        catch (...)
        {
            // Create a fresh exception
            reportException(new FATALERROR("Unhandled exception (not of type FatalError) in a parallel thread"));
        }
    }
}

////////////////////////////////////////////////////////////////////

void MultiProcessParallel::waitForChildThread()
{
    // Wait until the child thread is inactive
    std::unique_lock<std::mutex> lock(_mutex);
    while (_active) _conditionMain.wait(lock);
}

////////////////////////////////////////////////////////////////////

void MultiProcessParallel::reportException(FatalError* exception)
{
    // Need to lock, in case multiple threads throw simultaneously
    std::unique_lock<std::mutex> lock(_mutex);
    if (!_exception)  // only store the first exception thrown
    {
        _exception = exception;

        // Make the parent thread stop by taking away its work
        _maxIndex = 0;  // another thread will see either the old value, or zero
    }
}

///////////////////////////////////////////////////////////////////
