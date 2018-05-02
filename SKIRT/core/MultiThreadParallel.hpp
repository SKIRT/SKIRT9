/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef MULTITHREADPARALLEL_HPP
#define MULTITHREADPARALLEL_HPP

#include "Parallel.hpp"
#include <atomic>
#include <condition_variable>
#include <mutex>
#include <thread>
class FatalError;

////////////////////////////////////////////////////////////////////

/** This class implements the Parallel base class interface using multiple execution threads in a
    single process.

    When an exception is thrown by one of the threads excuting the target function, all other
    threads are gracefully shut down and a FatalError exception is thrown in the context of the
    thread that invoked the call() function. If the original exception was a FatalError instance,
    the newly thrown exception is a copy thereof. Otherwise a fresh FatalError instance is created
    with a generic error message.

    This class uses the standard low-level C++ multi-threading capabilities. It is designed to
    minimize the run-time overhead for handing out parallel tasks. Between invocations of the
    call() function, the parallel threads are put in wait so that they consume no CPU cycles (and
    very little memory).
*/
class MultiThreadParallel : public Parallel
{
    friend class ParallelFactory;       // so ParallelFactory can access our private constructor

    //============= Construction - Destruction =============

private:
    /** Constructs a MultiThreadParallel instance with the specified number of execution threads.
        The constructor is private; use the ParallelFactory::parallel() function instead. */
    explicit MultiThreadParallel(int threadCount);

public:
    /** Destructs the instance and its parallel threads. */
    ~MultiThreadParallel();

    //======================== Other Functions =======================

public:
    /** This function implements the call() interface described in the Parallel base class for the
        parallelization scheme offered by this subclass. */
    void call(std::function<void(size_t firstIndex, size_t numIndices)> target, size_t maxIndex) override;

private:
    /** The function that gets executed inside each of the parallel threads. */
    void run(int threadIndex);

    /** The function to do the actual work; used by call() and run(). */
    void doWork();

    /** A function to report an exception; used by doWork(). */
    void reportException(FatalError* exception);

    /** A function to wait for the parallel threads; used by constructor and call(). */
    void waitForThreads();

    /** This helper function returns true if at least one of the parallel threads (not including
        the parent thread) is still active, and false if not. This function does not perform any
        locking; callers should lock the shared data members of this class instance. */
    bool threadsActive();

    //======================== Data Members ========================

private:
    // data members keeping track of the threads
    int _threadCount;                   // the total number of threads, including the parent thread
    std::thread::id _parentThread;      // the ID of the thread that invoked our constructor
    std::vector<std::thread> _threads;  // the parallel threads (other than the parent thread)

    // synchronization
    std::mutex _mutex;                         // the mutex to synchronize with the parallel threads
    std::condition_variable _conditionExtra;   // the wait condition used by the parallel threads
    std::condition_variable _conditionMain;    // the wait condition used by the main thread

    // data members shared by all threads; changes are protected by a mutex
    std::function<void(size_t,size_t)> _target;    // the target function to be called
    size_t _chunkSize{0};               // the number of indices in all but the last chunk
    size_t _maxIndex{0};                // the maximum index (i.e. limiting the last chunk)
    FatalError* _exception{nullptr};    // a pointer to a heap-allocated copy of the exception thrown by a work thread
                                        // ... or zero if no exception was thrown
    std::vector<bool> _active;          // flag for each parallel thread (other than the parent thread)
                                        // ... that indicates whether the thread is currently active
    bool _terminate{false};             // becomes true when the parallel threads must exit

    // data member shared by all threads; incrementing is atomic (no need for protection)
    std::atomic<size_t> _nextIndex{0};  // tthe first index of the next availabe chunk
};

////////////////////////////////////////////////////////////////////

#endif
