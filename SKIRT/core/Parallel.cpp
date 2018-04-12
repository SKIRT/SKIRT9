/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "Parallel.hpp"
#include "FatalError.hpp"
#include "ParallelFactory.hpp"
#include "ProcessAssigner.hpp"

////////////////////////////////////////////////////////////////////

Parallel::Parallel(int threadCount, ParallelFactory* factory)
{
    // Cache the number of threads
    _threadCount = threadCount;

    // Remember the ID of the current thread
    _parentThread = std::this_thread::get_id();
    factory->addThreadIndex(_parentThread, 0);

    // Initialize shared data members and launch threads in a critical section
    {
        std::unique_lock<std::mutex> lock(_mutex);

        // Initialize shared data members
        _target = nullptr;
        _assigner = nullptr;
        _limit = 0;
        _active.assign(threadCount, true);
        _exception = nullptr;
        _terminate = false;
        _next = 0;

        // Create the extra parallel threads with one-based index (parent thread has index zero)
        for (int index = 1; index < threadCount; index++)
        {
            _threads.push_back(std::thread(&Parallel::run, this, index));
            factory->addThreadIndex(_threads.back().get_id(), index);
        }
    }

    // Wait until all parallel threads are ready
    waitForThreads();
}

////////////////////////////////////////////////////////////////////

Parallel::~Parallel()
{
    // Ask the parallel threads to exit in a critical section
    {
        std::unique_lock<std::mutex> lock(_mutex);

        // Ask the parallel threads to exit
        _terminate = true;
        _conditionExtra.notify_all();
    }

    // Wait for them to do so
    for (auto& thread: _threads) thread.join();
}

////////////////////////////////////////////////////////////////////

int Parallel::threadCount() const
{
    return _threadCount;
}

////////////////////////////////////////////////////////////////////

void Parallel::call(ParallelTarget* target, const ProcessAssigner* assigner, size_t repetitions)
{
    size_t assigned = assigner->assigned();
    call(target, assigner, assigned*repetitions, assigned);
}

////////////////////////////////////////////////////////////////////

void Parallel::call(ParallelTarget* target, size_t maxIndex, size_t repetitions)
{
    call(target, nullptr, maxIndex*repetitions, maxIndex);
}

////////////////////////////////////////////////////////////////////

void Parallel::call(ParallelTarget* target, const ProcessAssigner* assigner, size_t limit, size_t loopRange)
{
    // Verify that we're being called from our parent thread
    if (std::this_thread::get_id() != _parentThread)
        throw FATALERROR("Parallel call not invoked from thread that constructed this object");

    // Initialize shared data members and activate threads in a critical section
    {
        std::unique_lock<std::mutex> lock(_mutex);

        // Copy the arguments so they can be used from any of the threads
        _target = target;
        _assigner = assigner;
        _limit = limit;
        _loopRange = loopRange;

        // Initialize the number of active threads (i.e. not waiting for new work)
        _active.assign(_threadCount, true);

        // Clear the exception pointer
        _exception = 0;

        // Initialize the loop variable
        _next = 0;

        // Wake all parallel threads, if multithreading is allowed
        _conditionExtra.notify_all();
    }

    // Do some work ourselves as well
    doWork();

    // Wait until all parallel threads are done
    waitForThreads();

    // Check for and process the exception, if any
    if (_exception)
    {
        throw *_exception;  // throw by value (the memory for the heap-allocated exception is leaked)
    }
}

////////////////////////////////////////////////////////////////////

void Parallel::run(int threadIndex)
{
    while (true)
    {
        // Wait for new work in a critical section
        {
            std::unique_lock<std::mutex> lock(_mutex);

            // Indicate that this thread is no longer doing work
            _active[threadIndex] = false;

            // Tell the main thread when all parallel threads are inactive
            if (!threadsActive()) _conditionMain.notify_all();

            // Wait for new work
            while (true)
            {
                _conditionExtra.wait(lock);

                // Check for termination request (don't bother with _active; it's no longer used)
                if (_terminate) return;

                // Check that we actually have new work
                if (_active[threadIndex]) break;
            }
        }

        // Do work as long as some is available
        doWork();
    }
}

////////////////////////////////////////////////////////////////////

void Parallel::doWork()
{
    try
    {
        // Do work as long as some is available
        while (true)
        {
            size_t index = _next++;                  // get the next index atomically
            if (index >= _limit) break;              // break if no more are available

            // Repeat the same index range if necessary
            index = index % _loopRange;

            // Convert the index if using an assigner
            if (_assigner) index = _assigner->absoluteIndex(index);

            // Execute the body
            _target->body(index);
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

////////////////////////////////////////////////////////////////////

void Parallel::reportException(FatalError* exception)
{
    // Need to lock, in case multiple threads throw simultaneously
    std::unique_lock<std::mutex> lock(_mutex);
    if (!_exception)  // only store the first exception thrown
    {
        _exception = exception;

        // Make the other threads stop by taking away their work
        _limit = 0;  // this is safe because another thread will see either the old value, or zero
    }
}

////////////////////////////////////////////////////////////////////

void Parallel::waitForThreads()
{
    // Wait until all parallel threads are inactive
    std::unique_lock<std::mutex> lock(_mutex);
    while (threadsActive()) _conditionMain.wait(lock);
}

////////////////////////////////////////////////////////////////////

bool Parallel::threadsActive()
{
    // Check for active threads, skipping the parent thread with index 0
    return std::any_of(_active.begin()+1, _active.end(), [](bool flag){return flag;});
}

////////////////////////////////////////////////////////////////////
