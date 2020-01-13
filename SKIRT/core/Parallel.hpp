/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef PARALLEL_HPP
#define PARALLEL_HPP

#include "Basics.hpp"
#include <functional>

////////////////////////////////////////////////////////////////////

/** Parallel is an abstract base class for subclasses that implement various parallelization
    schemes using one or more execution threads and/or one or more processes in an MPI group.

    A Parallel subclass instance can be created only through the parallel() function of the
    ParallelFactory class. This function constructs the appropriate Parallel subclass instance
    depending on the requested task allocation mode and the available parallel resources. The
    client accesses the returned Parallel subclass instance using the common interface defined in
    this abstract base class.

    The call() function offered by this interface executes a specified target function \f$N\f$
    times as if it were part of a for loop over a range of indices from zero to \f$N-1\f$. Each
    index in the range represents a particular task. To reduce the overhead of handing out the
    tasks, the loop is actualy chopped into \em chunks of consecutive indices. Rather than a single
    index, the target function is handed the first index of the chunk and the number of indices
    (tasks) in the chunk, and it is expected to iterate over the specified index range. The chunk
    sizes are determined automatically to achieve optimal load balancing given the available
    parallel resources, while still maximally reducing the overhead of handing out the chunks.

    Because it may be invoked from multiple parallel threads, the target function should protect
    (write) access to shared data through the use of some synchronization mechanism. Shared data
    includes essentially everything other than local variables. Also note that this requirement
    also holds for all functions directly or indirectly called by the target function.

    A particular Parallel subclass instance can be reused for calling various target functions,
    reducing the overhead of creating and destroying the resources employed by the parallelization
    schem. One can also use multiple Parallel instances in a program. For example, a parallelized
    target function can invoke the call() function on another Parallel subclass instance that is
    requested from a ParallelFactory constructed and destructed within the scope of the target
    function. Recursively invoking the call() function on the same Parallel instance is not allowed
    and results in undefined behavior. */
class Parallel
{
    //============= Construction - Destruction =============

protected:
    /** The constructor is declared protected because this is an abstract class. At the level of
        this abstract class, the constructor does not do anything. */
    Parallel() {}

public:
    /** Destructs the instance and releases all resources. At the level of this abstract class, the
        destructor does not do anything. */
    virtual ~Parallel() {}

    //======================== Other Functions =======================

public:
    /** This function calls the specified target function repeatedly for index chunks that, taken
         together, will exactly cover the index range from zero to \em maxIndex-1. Each index chunk
         is specified to the target function through its two arguments, \em firstIndex and \em
         numIndices. The invocations of the target function will be distributed over the parallel
         resources in an unpredicable manner, and the various index chunks may be processed in
         arbitrary order.

         The index chunk sizes are determined automatically to achieve optimal load balancing given
         the available parallel resources, while still maximally reducing the overhead of handing
         out the chunks. */
    virtual void call(size_t maxIndex, std::function<void(size_t firstIndex, size_t numIndices)> target) = 0;
};

////////////////////////////////////////////////////////////////////

#endif
