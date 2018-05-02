/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef MULTIPROCESSPARALLEL_HPP
#define MULTIPROCESSPARALLEL_HPP

#include "Parallel.hpp"
#include <atomic>
#include <condition_variable>
#include <mutex>
#include <thread>
class FatalError;

////////////////////////////////////////////////////////////////////

/** This class implements the Parallel base class interface using a single thread in each of
    multiple processes. TODO: document and implement. */
class MultiProcessParallel : public Parallel
{
    friend class ParallelFactory;       // so ParallelFactory can access our private constructor

    //============= Construction - Destruction =============

private:
    /** Constructs a MultiProcessParallel instance. The specified number of execution threads is
        ignored. The number of processes is retrieved from the ProcessManager. In the root process,
        a child thread is created (and put on hold) so that the parent thread can be use to
        communicate with the other processes. This constructor is private; use the
        ParallelFactory::parallel() function instead. */
    explicit MultiProcessParallel(int threadCount);

public:
    /** Destructs the instance and, in the root process, its parallel child thread. */
    ~MultiProcessParallel();

    //======================== Other Functions =======================

public:
    /** This function implements the call() interface described in the Parallel base class for the
        parallelization scheme offered by this subclass. */
    void call(std::function<void(size_t firstIndex, size_t numIndices)> target, size_t maxIndex) override;

private:
    /** The function that gets executed inside the parallel child thread. Used only in root process. */
    void run();

    /** A function to wait for the parallel child thread. Used only in root process. */
    void waitForChildThread();

    /** A function to report an exception. Used only in root process. */
    void reportException(FatalError* exception);

    //======================== Data Members ========================

private:
    // data members are used only in the root process

    // child thread and synchronization
    std::thread _childThread;                   // the parallel child thread
    std::mutex _mutex;                          // the mutex to synchronize with the child thread
    std::condition_variable _conditionChild;    // the wait condition used by the child thread
    std::condition_variable _conditionMain;     // the wait condition used by the parent thread

    // data members shared by both threads; changes are protected by a mutex
    std::function<void(size_t,size_t)> _target; // the target function to be called
    size_t _chunkSize{0};                       // the number of indices in all but the last chunk
    size_t _maxIndex{0};                        // the maximum index (i.e. limiting the last chunk)
    FatalError* _exception{nullptr};            // a pointer to a heap-allocated copy of the exception thrown
                                                // ... by the child thread or zero if no exception was thrown
    bool _active{true};                         // flag that indicates whether the child thread is currently active
    bool _terminate{false};                     // becomes true when the child thread must exit

    // data member shared by both threads; incrementing is atomic (no need for protection)
    std::atomic<size_t> _nextIndex{0};          // the first index of the next availabe chunk
};

////////////////////////////////////////////////////////////////////

#endif
