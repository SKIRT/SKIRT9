/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef PARALLELFACTORY_HPP
#define PARALLELFACTORY_HPP

#include "SimulationItem.hpp"
#include "Parallel.hpp"
#include <thread>
#include <unordered_map>

/** A ParallelFactory object serves as a factory for instances of Parallel subclasses, called its
    children. An important property of a factory object is the maximum number of parallel execution
    threads per process to be handed to its children. During construction of a ParallelFactory
    instance, the maximum thread count is set to the return value of the defaultThreadCount()
    function. The value can be changed after construction by calling the setMaxThreadCount()
    function, but it should not be changed after any children have been requested. When requesting
    a new child, the client can specify a more stringent limit on the number of threads, but the
    factory's limit is never exceeded. This allows a client to request a Parallel object with a
    smaller number of threads, usually for performance reasons.

    A factory object assumes ownership for all its children. If a child of the appropriate type
    (see below) and with the appropriate number of threads already exists, it will be handed out
    again. As a result, a particular Parallel instance may be reused several times, reducing the
    overhead of creating and destroying the threads. However, the children of a particular factory
    should \em never be used in parallel. Also, recursively invoking the call() function on a
    Parallel instance is not allowed and results in undefined behavior. The recommended use is to
    have a single ParallelFactory instance per simulation, and to use yet another ParallelFactory
    instance to run multiple simulations at the same time.

    ParallelFactory clients can request a Parallel instance for one of the three task allocation
    modes described in the table below.

    Task mode | Description
    ----------|------------
    Distributed | All threads in all processes perform the tasks in parallel
    Duplicated | Each process performs all tasks; the results should be identical on all processes
    RootOnly | All threads in the root process perform the tasks in parallel; the other processes ignore the tasks

    In support of these task modes, the Parallel class has several subclasses, each implementing
    a specific parallelization scheme as described in the table below.

    Shorthand | %Parallel subclass | Description
    ----------|-----------------|------------
    S | SerialParallel | Single thread in the current process; isolated from any other processes
    MT | MultiThreadParallel | Multiple coordinated threads in the current process; isolated from any other processes
    MP | MultiProcessParallel | Single thread in each of multiple, coordinated processes
    MTP | HybridParallel | Multiple threads in each of multiple processes, all coordinated as a group
    0 | NullParallel | No operation; any requests for performing tasks are ignored

    Depending on the requested task mode and the current run-time configuration (number of processes
    and number of threads in each process), a ParallelFactory object hands out the appropriate
    Parallel subclass as listed in the table below. (Header: P=process, T=thread, 1=one, M=multiple.)

    Mode/Runtime | 1P 1T | 1P MT | MP 1T | MP MT |
    -------------|-------|-------|-------|-------|
    Distributed  |  S    |  MT   |  MP#  |  MTP# |
    Duplicated   |  S    |  MT   |  S    |  S*   |
    RootOnly     |  S    |  MT   |  S/0  |  MT/0 |

    (#) In Distributed mode with multiple processes, all threads require a different random number
        sequences. Therefore, the MultiProcessParallel and HybridParallel classes swith the Random
        instance associated with the simulation to abitrary mode before performing tasks, and back to
        predictable mode after performing the tasks.
    (*) In Duplicated mode with multiple processes, all tasks are performed by a single thread
        (in each process) because parallel threads executing tasks in an unpredictable order would
        see different random number sequences, possibly causing differences in the calculated results.
*/
class ParallelFactory : public SimulationItem
{
    friend class Parallel;

    //============= Construction - Setup - Destruction =============

public:
    /** This constructor creates a standalone parallel factory that is not hooked into a simulation
        item hierarchy; the caller is responsible for its destruction. Although the setup()
        function is not called by this constructor, no further setup is required to use the
        standalone object. */
    ParallelFactory();

    /** This constructor creates a parallel factory that is hooked up as a child to the specified
        parent in the simulation hierarchy, so that it will automatically be deleted. The setup()
        function is \em not called by this constructor. */
    explicit ParallelFactory(SimulationItem* parent);

    //====================== Other Functions =======================

public:
    /** Sets the maximum number of threads to be handed out to Parallel objects manufactured by
        this factory object. The minimum value is 1 thread. */
    void setMaxThreadCount(int value);

    /** Returns the maximum number of threads to be handed out to Parallel objects manufactured by
        this factory object. */
    int maxThreadCount() const;

    /** Returns the number of logical cores detected on the computer running the code, with a
        minimum of one and a maximum of 24 (additional threads in single process do not increase
        performance). */
    static int defaultThreadCount();

    /** Returns a Parallel instance with a particular number of execution threads. If the argument
        is zero or omitted, the number of threads equals the factory maximum. If the argument is
        nonzero, the number of threads is the smaller of the factory maximum and the specified
        maximum. */
    Parallel* parallel(int maxThreadCount=0);

    /** Returns an index for the parallel thread from which this function is called. When invoked
        from within a loop body being iterated by one of the factory's Parallel children, the
        function returns an index from zero to the number of threads in the Parallel instance minus
        one. When invoked from a thread that does not belong to any of the factory's children, the
        function throws a fatal error. */
    int currentThreadIndex() const;

private:
    /** Adds a dictionary item linking the specified thread to a particular index. This is a
        private function used from the Parallel() constructor to provide the information required
        by the currentThreadIndex() function. */
    void addThreadIndex(std::thread::id threadid, int index);

    //======================== Data Members ========================

private:
    // The maximum thread count for the factory, initialized to the default maximum number of threads
    int _maxThreadCount{ defaultThreadCount() };

    // The thread that invoked our constructor, initialized - obviously - upon construction
    std::thread::id _parentThread{ std::this_thread::get_id() };

    // The set of our children, keyed on number of threads; initially empty
    std::unordered_map<int, std::unique_ptr<Parallel>> _children;

    // The index for each child thread and for the parent thread; the latter is added upon construction
    std::unordered_map<std::thread::id, int> _indices{ { std::this_thread::get_id(), 0 } };
};

#endif
