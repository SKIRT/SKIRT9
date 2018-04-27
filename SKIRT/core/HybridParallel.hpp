/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef HYBRIDPARALLEL_HPP
#define HYBRIDPARALLEL_HPP

#include "Parallel.hpp"

////////////////////////////////////////////////////////////////////

/** This class implements the Parallel base class interface using multiple threads in each of
    multiple processes. TODO: document and implement. */
class HybridParallel : public Parallel
{
    friend class ParallelFactory;       // so ParallelFactory can access our private constructor

    //============= Construction - Destruction =============

private:
    /** Constructs a HybridParallel instance using the specified number of execution threads. The
        number of processes is retrieved from the ProcessManager. This constructor is private; use
        the ParallelFactory::parallel() function instead. */
    explicit HybridParallel(int threadCount);

    //======================== Other Functions =======================

public:
    /** This function implements the call() interface described in the Parallel base class for the
        parallelization scheme offered by this subclass. */
    void call(std::function<void(size_t firstIndex, size_t numIndices)> target, size_t maxIndex) override;
};

////////////////////////////////////////////////////////////////////

#endif
