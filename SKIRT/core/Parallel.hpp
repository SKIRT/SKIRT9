/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef PARALLEL_HPP
#define PARALLEL_HPP

#include "Basics.hpp"
#include <atomic>
#include <condition_variable>
#include <functional>
#include <mutex>
#include <thread>
class FatalError;
class ParallelFactory;

////////////////////////////////////////////////////////////////////

/** This class supports a simple parallel execution model similar to a for loop. In essence, the
    body of the for loop consists of a function specified by the caller which gets passed the index
    of the current iteration.

    To reduce the overhead of handing out the work and invoking the target, the loop is actualy
    carved into \em chunks of consecutive indices. Rather than a single index, the target function
    is handed the first index of the chunk and the number of indices in the chunk, and it is
    expected to iterate over the specified index range. By default, the chunk sizes are determined
    automatically to achieve optimal load balancing given the number of available parallel threads,
    while still maximally reducing the overhead of handing out the chunks. If desired, the user can
    override this behavior so that all chunks have a size of one, and the target function can
    ignore its second argument and does not need to iterate over the chunk size.

    A Parallel instance can be created only through the ParallelFactory class. The default
    construction determines a reasonable number of threads for the computer on which the code is
    running, and the call() function distributes the work over these threads.

    The parallelized body should protect (write) access to shared data through the use of some
    synchronization mechanism. Shared data includes the data members of the target object (in the
    context of which the target function executes) and essentially everything other than local
    variables. Also note that the parallelized body includes all functions directly or indirectly
    called by the target function.

    When an exception is thrown by one of the threads excuting the parallelized body, all other
    threads are gracefully shut down and a FatalError exception is thrown in the context of the
    thread that invoked the call() function. If the original exception was a FatalError instance,
    the newly thrown exception is a copy thereof. Otherwise a fresh FatalError instance is created
    with a generic error message.

    Between invocations of the call() function, the parallel threads are put in wait so that they
    consume no CPU cycles (and very little memory). Thus a particular Parallel instance can be
    reused many times for calling various member functions in various objects, reducing the
    overhead of creating and destroying the threads. One can also use multiple Parallel instances
    in an application. For example, a parallelized body can invoke the call() function on another
    Parallel instance that is constructed and destructed within the scope of the loop body.
    Recursively invoking the call() function on the same Parallel instance is not allowed and
    results in undefined behavior.

    The Parallel class uses C++ multi-threading capabilities, avoiding the complexities of using
    yet another library (such as OpenMP) and making it possible to properly handle exceptions. It
    is designed to minimize the run-time overhead for loops with many iterations. */
class Parallel
{
    friend class ParallelFactory;

    //============= Construction - Setup - Destruction =============

private:
    /** Constructs a Parallel instance with the specified number of execution threads. The
        constructor is not public; use the ParallelFactory::parallel() function instead. */
    Parallel(int threadCount, ParallelFactory* factory);

public:
    /** Destructs the instance and its parallel threads. */
    ~Parallel();

    //======================== Other Functions =======================

public:
    /** Returns the number of threads used by this instance. */
    int threadCount() const;

    /** Calls the specified target function repeatedly for index chunks that, taken together, will
        exactly cover the range from zero to \em maxIndex-1. Each index chunk is specified to the
        target function through its two arguments, \em firstIndex and \em numIndices. The
        invocations of the target function will be distributed over the parallel threads in an
        unpredicable manner, and the various chunks may be processed in arbitrary order.

        By default, the index chunk sizes are determined automatically to achieve optimal load
        balancing given the number of available parallel threads, while still maximally reducing
        the overhead of handing out the chunks. If the optional \em useChunksOfSizeOne flag is set
        to true, however, all chunks wil have a size of one. In this case, the target function can
        ignore its second argument and does not need to iterate over the chunk size. */
    void call(std::function<void(size_t firstIndex, size_t numIndices)> target,
              size_t maxIndex, bool useChunksOfSizeOne=false);

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
    size_t _numChunks;          // the number of chunks, or invocations of the target function
    size_t _chunkSize;          // the number of indices in all but the last chunk
    size_t _maxIndex;           // the maximum index (i.e. limiting the last chunk)
    std::vector<bool> _active;  // flag for each parallel thread (other than the parent thread)
                                // ... that indicates whether the thread is currently active
    FatalError* _exception;     // a pointer to a heap-allocated copy of the exception thrown by a work thread
                                // ... or zero if no exception was thrown
    bool _terminate;            // becomes true when the parallel threads must exit

    // data member shared by all threads; incrementing is atomic (no need for protection)
    std::atomic<size_t> _next;  // the current index of the for loop being implemented
};

////////////////////////////////////////////////////////////////////

#endif
