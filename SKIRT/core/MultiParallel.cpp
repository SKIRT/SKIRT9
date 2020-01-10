/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "MultiParallel.hpp"
#include "FatalError.hpp"

////////////////////////////////////////////////////////////////////

void MultiParallel::constructThreads(int numThreads)
{
    // Remember the number of threads
    _numThreads = numThreads;

    // Launch the child threads in a critical section
    {
        std::unique_lock<std::mutex> lock(_mutex);
        _active.assign(_numThreads, true);
        for (int index = 0; index != _numThreads; ++index)
        {
            _threads.push_back(std::thread(&MultiParallel::run, this, index));
        }
    }

    // Wait until all child threads are ready
    waitForThreads();
}

////////////////////////////////////////////////////////////////////

void MultiParallel::destroyThreads()
{
    // Ask the child threads to exit in a critical section
    {
        std::unique_lock<std::mutex> lock(_mutex);
        _terminate = true;
        _conditionChildren.notify_all();
    }

    // Wait for them to do so
    for (auto& thread : _threads) thread.join();
}

////////////////////////////////////////////////////////////////////

void MultiParallel::activateThreads()
{
    // Initialize shared data members and activate threads in a critical section
    std::unique_lock<std::mutex> lock(_mutex);
    _active.assign(_numThreads, true);
    _exception = nullptr;
    _conditionChildren.notify_all();
}

////////////////////////////////////////////////////////////////////

void MultiParallel::waitForThreads()
{
    // Wait until all parallel threads are inactive
    {
        std::unique_lock<std::mutex> lock(_mutex);
        while (threadsActive()) _conditionParent.wait(lock);
    }

    // Check for and process the exception, if any
    if (_exception)
    {
        throw *_exception;  // throw by value (the memory for the heap-allocated exception is leaked)
    }
}

////////////////////////////////////////////////////////////////////

void MultiParallel::run(int threadIndex)
{
    while (true)
    {
        // Wait for new work in a critical section
        {
            std::unique_lock<std::mutex> lock(_mutex);

            // Indicate that this thread is no longer doing work
            _active[threadIndex] = false;

            // Tell the main thread when all parallel threads are inactive
            if (!threadsActive()) _conditionParent.notify_all();

            // Wait for new work
            while (true)
            {
                _conditionChildren.wait(lock);

                // Check for termination request (don't bother with _active; it's no longer used)
                if (_terminate) return;

                // Check that we actually have new work
                if (_active[threadIndex]) break;
            }
        }

        // Do work as long as some is available for this cycle, and handle exceptions
        try
        {
            while (!_terminate && doSomeWork())
                ;
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

bool MultiParallel::threadsActive()
{
    // Check for active threads
    return std::any_of(_active.begin(), _active.end(), [](bool flag) { return flag; });
}

////////////////////////////////////////////////////////////////////

void MultiParallel::reportException(FatalError* exception)
{
    // Need to lock, in case multiple threads throw simultaneously
    std::unique_lock<std::mutex> lock(_mutex);
    if (!_exception)  // only store the first exception thrown
    {
        _exception = exception;

        // Make the other threads stop
        _terminate = true;
    }
}

////////////////////////////////////////////////////////////////////
