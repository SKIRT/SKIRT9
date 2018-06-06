/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "ProcessManager.hpp"
#include "FatalError.hpp"
#include <array>

#ifdef BUILD_WITH_MPI
#include <mpi.h>
#endif

////////////////////////////////////////////////////////////////////

int ProcessManager::_size{1};    // the number of processes: initialize to non-MPI default value
int ProcessManager::_rank{0};    // the rank of this process: initialize to non-MPI default value

////////////////////////////////////////////////////////////////////

#ifdef BUILD_WITH_MPI
namespace
{
    // Large messages will be broken up into pieces of the following size
    // (slightly under 2GB when data type is double)
    // because some MPI implementations dislike larger messages
    const size_t maxMessageSize = 250*1000*1000;
}
#endif

//////////////////////////////////////////////////////////////////////

void ProcessManager::initialize(int *argc, char ***argv)
{
#ifdef BUILD_WITH_MPI
    int initialized;
    MPI_Initialized(&initialized);
    if (!initialized)
    {
        // initialize MPI and verify that the implementation supports running multiple threads, as long
        // as we're calling MPI only from the main thread; this should avoid busy waits when blocking
        int provided = 0;
        MPI_Init_thread(argc, argv, MPI_THREAD_FUNNELED, &provided);
        if (provided < MPI_THREAD_FUNNELED)
            throw FATALERROR("MPI implementation does not support funneled threads");

        // get the process group size and our rank
        MPI_Comm_size(MPI_COMM_WORLD, &_size);
        MPI_Comm_rank(MPI_COMM_WORLD, &_rank);
    }
#else
    // the size and rank are statically initialized to the appropriate values
    (void)argc; (void)argv;
#endif
}

//////////////////////////////////////////////////////////////////////

void ProcessManager::finalize()
{
#ifdef BUILD_WITH_MPI
    MPI_Finalize();
#endif
}

void ProcessManager::abort(int exitcode)
{
#ifdef BUILD_WITH_MPI
    if (isMultiProc()) MPI_Abort(MPI_COMM_WORLD, exitcode);
#else
    (void)exitcode;
#endif
}

//////////////////////////////////////////////////////////////////////

namespace
{
    void throwInvalidChunkInvocation()
    {
        throw FATALERROR("Multi-process chunk allocation function called from inappropriate process");
    }
}

//////////////////////////////////////////////////////////////////////

bool ProcessManager::requestChunk(size_t& firstIndex, size_t& numIndices)
{
#ifdef BUILD_WITH_MPI
    if (isRoot()) throwInvalidChunkInvocation();

    std::array<int,1> sendbuf{{_rank}};  // we pass our rank so that the receiver can ignore MPI status
    std::array<size_t,2> recvbuf{{0,0}};
    MPI_Sendrecv(sendbuf.begin(), sendbuf.size(), MPI_INT, 0, 1,
                 recvbuf.begin(), recvbuf.size(), MPI_UNSIGNED_LONG, 0, 1,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    firstIndex = recvbuf[0];
    numIndices = recvbuf[1];
    return numIndices > 0;
#else
    (void)firstIndex; (void)numIndices;
    throwInvalidChunkInvocation();
    return false;
#endif
}

//////////////////////////////////////////////////////////////////////

int ProcessManager::waitForChunkRequest()
{
#ifdef BUILD_WITH_MPI
    if (!isMultiProc() || !isRoot()) throwInvalidChunkInvocation();

    std::array<int,1> recvbuf{{0}};
    MPI_Recv(recvbuf.begin(), recvbuf.size(), MPI_INT, MPI_ANY_SOURCE, 1,
             MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    return recvbuf[0];
#else
    throwInvalidChunkInvocation();
    return 0;
#endif
}

//////////////////////////////////////////////////////////////////////

void ProcessManager::serveChunkRequest(int rank, size_t firstIndex, size_t numIndices)
{
#ifdef BUILD_WITH_MPI
    if (!isMultiProc() || !isRoot()) throwInvalidChunkInvocation();

    std::array<size_t,2> sendbuf{{firstIndex,numIndices}};
    MPI_Send(sendbuf.begin(), sendbuf.size(), MPI_UNSIGNED_LONG, rank, 1,
             MPI_COMM_WORLD);
#else
    (void)rank; (void)firstIndex; (void)numIndices;
    throwInvalidChunkInvocation();
#endif
}

//////////////////////////////////////////////////////////////////////

void ProcessManager::wait()
{
#ifdef BUILD_WITH_MPI
    if (isMultiProc()) MPI_Barrier(MPI_COMM_WORLD);
#endif
}

//////////////////////////////////////////////////////////////////////

void ProcessManager::sumToAll(Array& arr)
{
#ifdef BUILD_WITH_MPI
    if (isMultiProc())
    {
        double* data = begin(arr);
        size_t remaining = arr.size();
        while (remaining > maxMessageSize)
        {
            MPI_Allreduce(MPI_IN_PLACE, data, maxMessageSize, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            data += maxMessageSize;
            remaining -= maxMessageSize;
        }
        if (remaining)
        {
            MPI_Allreduce(MPI_IN_PLACE, data, remaining, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        }
    }
#else
    (void)arr;
#endif
}

//////////////////////////////////////////////////////////////////////

void ProcessManager::sumToRoot(Array& arr)
{
#ifdef BUILD_WITH_MPI
    if (isMultiProc())
    {
        double* data = begin(arr);
        size_t remaining = arr.size();

        while (remaining > maxMessageSize)
        {
            if (isRoot())
                MPI_Reduce(MPI_IN_PLACE, data, maxMessageSize, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
            else
                MPI_Reduce(data, data, maxMessageSize, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

            remaining -= maxMessageSize;
            data += maxMessageSize;
        }
        if (remaining)
        {
            if (isRoot())
                MPI_Reduce(MPI_IN_PLACE, data, remaining, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
            else
                MPI_Reduce(data, data, remaining, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        }
    }
#else
    (void)arr;
#endif
}

//////////////////////////////////////////////////////////////////////
