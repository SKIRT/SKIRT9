/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef PARALLEL_HPP
#define PARALLEL_HPP

#include "ParallelTarget.hpp"
#include <atomic>
#include <condition_variable>
#include <mutex>
#include <thread>
class FatalError;
class ParallelFactory;
class ProcessAssigner;

////////////////////////////////////////////////////////////////////

/** This class supports a simple parallel execution model similar to a for loop. The body of the
    for loop consists of the body() function of a ParallelTarget subclass, or alternatively of
    some class member function (specified by the caller) with a single integer argument, which
    gets passed the index of the current iteration.

    For example, to invoke a member function 1000 times one would write something like:
        \code
        void Test::body(size_t index) { ... }
        ...
        SequentialAssigner assigner;
        ParallelFactory factory;
        ...
        assigner.assign(1000);
        factory.parallel()->call(this, &Test::body, &assigner);
        \endcode

    A Parallel instance can be created only through the ParallelFactory class. The default
    construction shown above determines a reasonable number of threads for the computer on which
    the code is running, and the call() function distributes the work over these threads.

    The parallelized body should protect (write-) access to shared data through the use of an
    std::mutex object. Shared data includes the data members of the target object (in the context
    of which the target function executes) and essentially everything other than local variables.
    Also note that the parallelized body includes all functions directly or indirectly called by
    the target function.

    When an exception is thrown by one of the threads excuting the parallelized body, all other
    threads are gracefully shut down and a FatalError exception is thrown in the context of the
    thread that invoked the call() function. If the original exception was a FatalError instance, the
    newly thrown exception is a copy thereof. Otherwise a fresh FatalError instance is created with
    a generic error message.

    Between invocations of the call() function, the parallel threads are put in wait so that
    they consume no CPU cycles (and very little memory). Thus a
    particular Parallel instance can be reused many times for calling various member functions in
    various objects, reducing the overhead of creating and destroying the threads. One can also
    use multiple Parallel instances in an application. For example, a parallelized body can invoke
    the call() function on another Parallel instance that is constructed and destructed within the
    scope of the loop body. Recursively invoking the call() function on the same Parallel instance is
    not allowed and results in undefined behavior.

    The Parallel class uses C++11 multi-threading capabilities, avoiding the complexities of using
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

    /** Calls the body() function of the specified specified target object a certain number of
        times, with the \em index argument of that function taking values that are determined by
        the \em assigner, which is also passed to this function. While the values assigned to a
        particular process are fixed at the moment this function is called, the work will be
        distributed over the parallel threads in an unpredicable manner. If the \em repetitions
        argument is > 1, this function will loop through the indices specified by the assigner
        multiple times. */
    void call(ParallelTarget* target, const ProcessAssigner* assigner, size_t repetitions = 1);

    /** Calls the body() function of the specified specified target object a certain number of
        times, with the \em index argument of that function taking values from 0 to \em maxIndex-1.
        The work will be distributed over the parallel threads in an unpredicable manner. If the
        \em repetitions argument is > 1, the loop over the indices will be repeated multiple times.
        With every repetition, the index starts from 0 and goes up to \em maxIndex-1. */
    void call(ParallelTarget* target, size_t maxIndex, size_t repetitions = 1);

    /** Calls the specified member function for the specified target object a certain number of
        times, with the \em index argument of that function taking values that are determined by
        the \em assigner, which is also passed to this function. While the values assigned to a
        particular process are fixed at the moment this function is called, the work will be
        distributed over the parallel threads in an unpredicable manner. If the \em repetitions
        argument is > 1, this function will loop through the indices specified by the assigner
        multiple times. */
    template<class T> void call(T* targetObject, void (T::*targetMember)(size_t index),
                                const ProcessAssigner* assigner, size_t repetitions = 1);

    /** Calls the specified member function for the specified target object a certain number of
        times, with the \em index argument of that function taking values from 0 to \em maxIndex-1.
        The work will be distributed over the parallel threads in an unpredicable manner. If the
        \em repetitions argument is > 1, the loop over the indices will be repeated multiple times.
        With every repetition, the index starts from 0 and goes up to \em maxIndex-1. */
    template<class T> void call(T* targetObject, void (T::*targetMember)(size_t index),
                                size_t maxIndex, size_t repetitions = 1);

private:
    /** This function gets called by all other versions of the call() function. It sets the data
        members shared by the threads and starts the parallel execution. */
    void call(ParallelTarget* target, const ProcessAssigner* assigner, size_t limit, size_t loopRange);

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

    //======================== Nested Classes =======================

private:
    /** The declaration for this template class is nested in the Parallel class declaration. It is
        used in the implementation of the call() template function to allow specifying a loop body
        in the form of an arbitrary target member function and target object. An instance of this
        target-type-dependent class bundles the relevant information into a single object with a
        type-independent base class (ParallelTarget) so that it can be passed to the regular (i.e.
        non-template) call() function. */
    template<class T> class Target : public ParallelTarget
    {
    public:
        /** Constructs a ParallelTarget instance with a body() function that calls the specified
            target member function on the specified target object. */
        Target(T* targetObject, void (T::*targetMember)(size_t index))
            : _targetObject(targetObject), _targetMember(targetMember) { }

        /** Calls the target member function on the target object specified in the constructor. */
        void body(size_t index) { (_targetObject->*(_targetMember))(index); }

    private:
        T* _targetObject;
        void (T::*_targetMember)(size_t index);
    };

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
    ParallelTarget* _target;    // the target to be called
    const ProcessAssigner* _assigner; // the process assigner
    size_t _limit;              // the total number of calls to the given function
    size_t _loopRange;          // the number of function calls in one iteration of the for loop
    std::vector<bool> _active;  // flag for each parallel thread (other than the parent thread)
                                // ... that indicates whether the thread is currently active
    FatalError* _exception;     // a pointer to a heap-allocated copy of the exception thrown by a work thread
                                // ... or zero if no exception was thrown
    bool _terminate;            // becomes true when the parallel threads must exit

    // data member shared by all threads; incrementing is atomic (no need for protection)
    std::atomic<size_t> _next;   // the current index of the for loop being implemented
};

////////////////////////////////////////////////////////////////////

// outermost portion of the call() template function implementation
template<class T> void Parallel::call(T* targetObject, void (T::*targetMember)(size_t index),
                                      const ProcessAssigner* assigner, size_t repetitions)
{
    Target<T> target(targetObject, targetMember);
    call(&target, assigner, repetitions);
}

////////////////////////////////////////////////////////////////////

template<class T> void Parallel::call(T* targetObject, void (T::*targetMember)(size_t index),
                                      size_t maxIndex, size_t repetitions)
{
    Target<T> target(targetObject, targetMember);
    call(&target, maxIndex, repetitions);
}

////////////////////////////////////////////////////////////////////

#endif
