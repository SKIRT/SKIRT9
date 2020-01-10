/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef CHUNKMAKER_HPP
#define CHUNKMAKER_HPP

#include "Basics.hpp"
#include <atomic>
#include <functional>

//////////////////////////////////////////////////////////////////////

/** The ChunkMaker class implements a heuristic for chopping a range of indices from zero to
    \f$N-1\f$ into smaller ranges of consecutive indices called \em chunks. An instance of the
    ChunkMaker class is used, for example, by some Parallel subclasses to distribute a range of
    tasks over parallel execution threads and/or processes. Assuming that each index maps to a
    particular task, the function performing the tasks is repeatedly handed the first index of the
    chunk and the number of indices in the chunk, and it is expected to iterate over the specified
    index range. The chunk sizes are determined by the heuristic in the ChunkMaker object to
    achieve optimal load balancing given the available parallel resources, while still maximally
    reducing the overhead of handing out the chunks. */
class ChunkMaker
{
public:
    /** The default (and only) constructor initializes the ChunkMaker object to an empty range. */
    ChunkMaker();

    /** This function initializes the ChunkMaker object to the specified range (from zero to
        \f$N-1\f$), using the specified number of threads and processes to help determine an
        appropriate chunk size. */
    void initialize(size_t maxIndex, int numThreads, int numProcs = 1);

    /** This function gets the next chunk, in the form of the first index and the number of indices
        in the chunk. If a chunk is still available, the function places a chunk index range in its
        arguments and returns true. If no more chunks are available, the output arguments remain
        unchanged and the function returns false. This function uses an atomic operation to obtain
        the next chunk so it can safely be called from multiple concurrent execution threads. */
    bool next(size_t& firstIndex, size_t& numIndices);

    /** This function gets the next chunk, and if one is still available, it calls the specified
        target with the corresponding first index and number of indices, and returns true. If no
        more chunks are available, the target is not invoked and this function returns false. This
        function uses an atomic operation to obtain the next chunk so it can safely be called from
        multiple concurrent execution threads, as long as the target function is thread-safe as
        well. */
    bool callForNext(const std::function<void(size_t firstIndex, size_t numIndices)>& target);

private:
    size_t _chunkSize{0};               // the number of indices in all but the last chunk
    size_t _maxIndex{0};                // the maximum index (i.e. limiting the last chunk)
    std::atomic<size_t> _nextIndex{0};  // the first index of the next available chunk
};

//////////////////////////////////////////////////////////////////////

#endif
