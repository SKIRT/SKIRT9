/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef SERIALPARALLEL_HPP
#define SERIALPARALLEL_HPP

#include "Parallel.hpp"

////////////////////////////////////////////////////////////////////

/** This class implements the Parallel base class interface using a single execution threads in a
    single process. In other words, all tasks are serialized in the thread invoking the call()
    function. Because of the simplicity of the "parallelization" scheme, overhead is minimal. Also,
    any exceptions thrown by the target function are simply passed through. */
class SerialParallel : public Parallel
{
    friend class ParallelFactory;  // so ParallelFactory can access our private constructor

    //============= Construction - Destruction =============

private:
    /** Constructs a SerialParallel instance. The specified number of execution threads is ignored.
        The constructor is private; use the ParallelFactory::parallel() function instead. */
    explicit SerialParallel(int threadCount);

    //======================== Other Functions =======================

public:
    /** This function implements the call() interface described in the Parallel base class for the
        parallelization scheme offered by this subclass. */
    void call(size_t maxIndex, std::function<void(size_t firstIndex, size_t numIndices)> target) override;
};

////////////////////////////////////////////////////////////////////

#endif
