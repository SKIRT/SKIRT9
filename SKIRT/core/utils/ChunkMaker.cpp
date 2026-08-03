/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "ChunkMaker.hpp"

//////////////////////////////////////////////////////////////////////

ChunkMaker::ChunkMaker() {}

//////////////////////////////////////////////////////////////////////

void ChunkMaker::initialize(size_t maxIndex, int numThreads, int numProcs)
{
    // Determine the chunk size
    const size_t numChunksPerThread = 8;  // empirical multiplicator to achieve acceptable load balancing
    _chunkSize = max(static_cast<size_t>(1), maxIndex / (numThreads * numProcs * numChunksPerThread));

    // Initialize the other data members
    _maxIndex = maxIndex;
    _nextIndex = 0;
}

//////////////////////////////////////////////////////////////////////

bool ChunkMaker::next(size_t& firstIndex, size_t& numIndices)
{
    size_t first = _nextIndex.fetch_add(_chunkSize);
    if (first < _maxIndex)
    {
        firstIndex = first;
        numIndices = min(_chunkSize, _maxIndex - first);
        return true;
    }
    return false;
}

//////////////////////////////////////////////////////////////////////

bool ChunkMaker::callForNext(const std::function<void(size_t, size_t)>& target)
{
    size_t first = _nextIndex.fetch_add(_chunkSize);
    if (first < _maxIndex)
    {
        target(first, min(_chunkSize, _maxIndex - first));
        return true;
    }
    return false;
}

//////////////////////////////////////////////////////////////////////
