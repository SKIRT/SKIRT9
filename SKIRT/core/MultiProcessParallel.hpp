/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef MULTIPROCESSPARALLEL_HPP
#define MULTIPROCESSPARALLEL_HPP

#include "MultiParallel.hpp"
#include "ChunkMaker.hpp"

////////////////////////////////////////////////////////////////////

/** This class implements the Parallel base class interface using a single execution thread in each
    of multiple processes. In the root process, the class in fact employs a second thread to
    perform the actual work, so that the parent thread can serve chunk requests to the other
    processes (the MPI functions should be called only from the main thread). This extra thread is
    not counted towards the number of threads specified by the user because it does not consume
    significant resources. In the other (non-root) processes, there is no extra thread; the work is
    performed in a loop that simply requests and performs new chunks.

    This class uses the facilities offered by the MultiParallel base class. */
class MultiProcessParallel : public MultiParallel
{
    friend class ParallelFactory;       // so ParallelFactory can access our private constructor

    //============= Construction - Destruction =============

private:
    /** Constructs a MultiProcessParallel instance. The specified number of execution threads is
        ignored. The number of processes is retrieved from the ProcessManager. In the root process,
        a child thread is created (and put on hold) so that the parent thread can be used to
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
    void call(size_t maxIndex, std::function<void(size_t firstIndex, size_t numIndices)> target) override;

private:
    /** The function to do the actual work, one chunk at a time. */
    bool doSomeWork() override;

    //======================== Data Members ========================

private:
    // data members are used only in the root process
    std::function<void(size_t,size_t)> _target; // the target function to be called
    ChunkMaker _chunkMaker;                     // the chunk maker
};

////////////////////////////////////////////////////////////////////

#endif
