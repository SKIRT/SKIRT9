/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef MULTIPROCESSPARALLEL_HPP
#define MULTIPROCESSPARALLEL_HPP

#include "Parallel.hpp"

////////////////////////////////////////////////////////////////////

/** This class implements the Parallel base class interface using a single thread in each of
    multiple processes. TODO: document and implement. */
class MultiProcessParallel : public Parallel
{
    friend class ParallelFactory;       // so ParallelFactory can access our private constructor

    //============= Construction - Destruction =============

private:
    /** Constructs a MultiProcessParallel instance. The specified number of execution threads is
        ignored. The number of processes is retrieved from the ProcessManager. This constructor is
        private; use the ParallelFactory::parallel() function instead. */
    explicit MultiProcessParallel(int threadCount);

    //======================== Other Functions =======================

public:
    /** This function implements the call() interface described in the Parallel base class for the
        parallelization scheme offered by this subclass. */
    void call(std::function<void(size_t firstIndex, size_t numIndices)> target, size_t maxIndex) override;
};

////////////////////////////////////////////////////////////////////

#endif
