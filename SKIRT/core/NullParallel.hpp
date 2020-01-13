/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef NULLPARALLEL_HPP
#define NULLPARALLEL_HPP

#include "Parallel.hpp"

////////////////////////////////////////////////////////////////////

/** This class implements the Parallel base class interface as a "no operation". In other words, it
    ignores any requests for performing tasks. This no-operation capability can be used, for
    example, to skip a certain body of work on all processeses except for the root process. */
class NullParallel : public Parallel
{
    friend class ParallelFactory;  // so ParallelFactory can access our private constructor

    //============= Construction - Destruction =============

private:
    /** Constructs a NullParallel instance. The specified number of execution threads is ignored.
        The constructor is private; use the ParallelFactory::parallel() function instead. */
    explicit NullParallel(int threadCount);

    //======================== Other Functions =======================

public:
    /** This function implements the call() interface described in the Parallel base class for the
        parallelization scheme offered by this subclass. */
    void call(size_t maxIndex, std::function<void(size_t firstIndex, size_t numIndices)> target) override;
};

////////////////////////////////////////////////////////////////////

#endif
