/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef MULTIPARALLEL_HPP
#define MULTIPARALLEL_HPP

#include "Parallel.hpp"
#include <atomic>
#include <condition_variable>
#include <mutex>
#include <thread>
class FatalError;

////////////////////////////////////////////////////////////////////

/** MultiParallel is an intermediate, abstract class that offers facilities for implementing
    parallel subclasses that use multiple parallel execution threads.

    Specifically, the class offers functions to create a given number of child threads and hand out
    work to them. The parent thread is never used to perform actual work, so that it remains
    available for communicating with other processes. Because the parent thread is supposed to
    consume very little resources, it is not counted towards the number of threads configured by
    the user.

    When an exception is thrown by one of the child threads, all other threads are gracefully shut
    down and a FatalError exception is thrown in the context of the parent thread. If the original
    exception was a FatalError instance, the newly thrown exception is a copy thereof. Otherwise a
    fresh FatalError instance is created with a generic error message.

    This class uses the standard low-level C++ multi-threading capabilities. It is designed to
    minimize the run-time overhead for handing out parallel tasks. Between invocations of the
    call() function, the parallel threads are put in wait so that they consume no CPU cycles (and
    very little memory). */
class MultiParallel : public Parallel
{
    //============== Facilities offered by this class ==============

protected:
    /** This function constructs the specified number of parallel child threads (not including the
        parent thread) and waits for them to become ready (in the inactive state). */
    void constructThreads(int numThreads);

    /** This function destructs the child threads constucted with the constructThreads() function.
        */
    void destroyThreads();

    /** This function activates the child threads so they start doing work (i.e. calling the
        doSomeWork() function until it returns false). */
    void activateThreads();

    /** This function blocks until all child threads have become inactive (i.e. the doSomeWork()
        function has returned false for all threads). */
    void waitForThreads();

    /** This function returns the number of parallel child threads (not including the parent
        thread) specified to constructThreads(). */
    int numThreads() { return _numThreads; }

private:
    /** This function gets executed inside each of the parallel threads. */
    void run(int threadIndex);

    /** This function returns true if at least one of the child threads is still active, and false
        if not. This function does not perform any locking; callers should lock the shared data
        members of this class instance. */
    bool threadsActive();

    /** This function reports an exception from within any thread. */
    void reportException(FatalError* exception);

    //========= Functions to be implemented by subclasses ==========

    /** The function to do the actual work; called from within run(). The function should perform
        some limited amount of work and then return true if more work might be available for this
        cycle, and false if not. */
    virtual bool doSomeWork() = 0;

    //======================== Data Members ========================

private:
    // the threads
    int _numThreads{0};                 // the number of child threads (not including the parent thread)
    std::vector<std::thread> _threads;  // the child threads

    // synchronization
    std::mutex _mutex;                           // the mutex to synchronize the threads
    std::condition_variable _conditionChildren;  // the wait condition used by the child threads
    std::condition_variable _conditionParent;    // the wait condition used by the parent thread

    // data members shared by all threads; changes are protected by a mutex
    std::vector<bool> _active;            // flag for each child thread
                                          // ... that indicates whether the thread is currently active
    FatalError* _exception{nullptr};      // a pointer to a heap-allocated copy of the exception thrown
                                          // ...  by a child thread or null if no exception was thrown
    std::atomic<bool> _terminate{false};  // becomes true when the child threads must exit
};

////////////////////////////////////////////////////////////////////

#endif
