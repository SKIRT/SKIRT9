/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef PARALLELFACTORY_HPP
#define PARALLELFACTORY_HPP

#include "SimulationItem.hpp"
#include <map>
#include <thread>
class Parallel;

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

    ParallelFactory clients can request a Parallel instance for one of the two task allocation
    modes described in the table below.

    Task mode | Description
    ----------|------------
    Distributed | All threads in all processes perform the tasks in parallel
    RootOnly | All threads in the root process perform the tasks in parallel; the other processes ignore the tasks

    In support of these task modes, the Parallel class has several subclasses, each implementing
    a specific parallelization scheme as described in the table below.

    Shorthand | %Parallel subclass | Description
    ----------|-----------------|------------
    S | SerialParallel | Single thread in the current process; isolated from any other processes
    MT | MultiThreadParallel | Multiple coordinated threads in the current process; isolated from any other processes
    MTP | MultiHybridParallel | One or more threads in each of multiple processes, all coordinated as a group
    0 | NullParallel | No operation; any requests for performing tasks are ignored

    Depending on the requested task mode and the current run-time configuration (number of processes
    and number of threads in each process), a ParallelFactory object hands out the appropriate
    Parallel subclass as listed in the table below. (Header: P=process, T=thread, 1=one, M=multiple.)

    Mode/Runtime | 1P 1T | 1P MT | MP 1T | MP MT |
    -------------|-------|-------|-------|-------|
    Distributed  |  S    |  MT   |  MTP  |  MTP  |
    RootOnly     |  S    |  MT   |  S/0  |  MT/0 |

*/
class ParallelFactory : public SimulationItem
{
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

    /** The destructor releases any resources held by the parallel factory. */
    ~ParallelFactory();

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

    /** This enumeration includes a constant for each task allocation mode supported by ParallelFactory
     * and the Parallel subclasses. */
    enum class TaskMode { Distributed, RootOnly };

    /** This function returns a Parallel subclass instance of the appropriate type and with an
        appropriate number of execution threads, depending on the requested task allocation mode,
        the specified maximum number of threads, and the current run-time environment (number of
        threads in the factory maximum and number of processes in the MPI group). The recipe for
        determining the appropriate paralllization scheme is described in the main documentation
        for this class.

        The first argument specifies the task mode. The second argument, if present, limits the
        maximum number of threads to the specified number. */
    Parallel* parallel(TaskMode mode, int maxThreadCount = 0);

    /** This function calls the parallel() function for the Distributed task allocation mode. */
    Parallel* parallelDistributed(int maxThreadCount = 0) { return parallel(TaskMode::Distributed, maxThreadCount); }

    /** This function calls the parallel() function for the RootOnly task allocation mode. */
    Parallel* parallelRootOnly(int maxThreadCount = 0) { return parallel(TaskMode::RootOnly, maxThreadCount); }

    //======================== Data Members ========================

private:
    // The maximum thread count for the factory, initialized to the default maximum number of threads
    int _maxThreadCount{defaultThreadCount()};

    // The thread that invoked our constructor, initialized - obviously - upon construction
    std::thread::id _parentThread{std::this_thread::get_id()};

    // Private enumeration of the supported Parallel subclasses
    enum class ParallelType { Null = 0, Serial, MultiThread, MultiHybrid };

    // The collection of our children, keyed on Parallel subclass type and number of threads; initially empty
    std::map<std::pair<ParallelType, int>, std::unique_ptr<Parallel>> _children;
};

#endif
