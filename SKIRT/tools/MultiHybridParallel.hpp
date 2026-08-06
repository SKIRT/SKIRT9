/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef MULTIHYBRIDPARALLEL_HPP
#define MULTIHYBRIDPARALLEL_HPP

#include "ChunkMaker.hpp"
#include "MultiParallel.hpp"

////////////////////////////////////////////////////////////////////

/** This class implements the Parallel base class interface using multiple threads in each of
    multiple processes. In each process, the actual work is performed in child threads created for
    that purpose, while the parent thread is used for communication among the processes (the MPI
    functions should be called only from the main thread). In the root process, the parent thread
    serves chunks of work to the other processes. In the non-root processes, the parent thread
    requests chunks from the root process for all of the local parallel threads. The extra thread
    in each process is not counted towards the number of threads specified by the user because the
    communication does not consume significant resources.

    This class uses the facilities offered by the MultiParallel base class. */
class MultiHybridParallel : public MultiParallel
{
    friend class ParallelFactory;  // so ParallelFactory can access our private constructor

    //============= Construction - Destruction =============

private:
    /** Constructs a HybridParallel instance using the specified number of execution threads. The
        number of processes is retrieved from the ProcessManager. In each process, the specified
        number of child threads is created (and put on hold) so that the parent thread can be used
        to communicate with the other processes. This constructor is private; use the
        ParallelFactory::parallel() function instead. */
    explicit MultiHybridParallel(int threadCount);

public:
    /** Destructs the instance and its parallel child threads. */
    ~MultiHybridParallel();

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
    // used in all processes
    std::function<void(size_t, size_t)> _target;  // the target function to be called

    // used only in the root process; shared between threads
    ChunkMaker _chunkMaker;  // the chunk maker

    // used only in non-root processes; shared between threads
    std::mutex _mutex;                           // the mutex to synchronize the threads
    std::condition_variable _conditionChildren;  // the wait condition used by the child threads
    std::condition_variable _conditionParent;    // the wait condition used by the parent thread
    int _requests{0};                            // the number of outstanding chunk requests from child threads
    bool _ready{false};                          // true if firstIndex/numIndices represent a valid, unconsumed chunk
    bool _done{false};                           // true if there are no more chunks to be served
    size_t _firstIndex{0};                       // the first index of the new chunk being served
    size_t _numIndices{0};                       // the number of indices of the new chunk being served
};

////////////////////////////////////////////////////////////////////

#endif
