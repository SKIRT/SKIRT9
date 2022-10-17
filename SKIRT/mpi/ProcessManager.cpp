/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "ProcessManager.hpp"
#include "FatalError.hpp"
#include <array>

#ifdef BUILD_WITH_MPI
#    include <mpi.h>
#    include <chrono>
#    include <thread>
#endif

////////////////////////////////////////////////////////////////////

int ProcessManager::_size{1};  // the number of processes: initialize to non-MPI default value
int ProcessManager::_rank{0};  // the rank of this process: initialize to non-MPI default value

////////////////////////////////////////////////////////////////////

#ifdef BUILD_WITH_MPI
namespace
{
    // Large messages will be broken up into pieces of the following size
    // (slightly under 2GB when data type is double)
    // because some MPI implementations dislike larger messages
    const size_t maxMessageSize = 250 * 1000 * 1000;
}
#endif

//////////////////////////////////////////////////////////////////////

void ProcessManager::initialize(int* argc, char*** argv)
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
        if (provided < MPI_THREAD_FUNNELED) throw FATALERROR("MPI implementation does not support funneled threads");

        // get the process group size and our rank
        MPI_Comm_size(MPI_COMM_WORLD, &_size);
        MPI_Comm_rank(MPI_COMM_WORLD, &_rank);
    }
#else
    // the size and rank are statically initialized to the appropriate values
    (void)argc;
    (void)argv;
#endif
}

//////////////////////////////////////////////////////////////////////

void ProcessManager::finalize()
{
#ifdef BUILD_WITH_MPI
    MPI_Finalize();
#endif
}

//////////////////////////////////////////////////////////////////////

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
    std::function<void(string)> _logger;  // the usage logger; intialize to empty
}

//////////////////////////////////////////////////////////////////////

void ProcessManager::setLogger(std::function<void(string)> logger)
{
    _logger = logger;
}

//////////////////////////////////////////////////////////////////////

void ProcessManager::clearLogger()
{
    _logger = nullptr;
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

    if (_logger) _logger("MPI BEGIN: request chunk");
    std::array<int, 1> sendbuf{{_rank}};  // we pass our rank so that the receiver can ignore MPI status
    std::array<size_t, 2> recvbuf{{0, 0}};
    MPI_Sendrecv(sendbuf.begin(), sendbuf.size(), MPI_INT, 0, 1, recvbuf.begin(), recvbuf.size(), MPI_UNSIGNED_LONG, 0,
                 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    firstIndex = recvbuf[0];
    numIndices = recvbuf[1];
    if (_logger) _logger("MPI END: received chunk " + std::to_string(firstIndex) + ", " + std::to_string(numIndices));
    return numIndices > 0;
#else
    (void)firstIndex;
    (void)numIndices;
    throwInvalidChunkInvocation();
    return false;
#endif
}

//////////////////////////////////////////////////////////////////////

int ProcessManager::waitForChunkRequest()
{
#ifdef BUILD_WITH_MPI
    if (!isMultiProc() || !isRoot()) throwInvalidChunkInvocation();

    if (_logger) _logger("MPI BEGIN: wait for chunk request");
    while (true)  // avoid using CPU while waiting for a message
    {
        int flag;
        MPI_Iprobe(MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, &flag, MPI_STATUS_IGNORE);
        if (flag) break;
        std::this_thread::sleep_for(std::chrono::milliseconds(20));
    }
    std::array<int, 1> recvbuf{{0}};
    MPI_Recv(recvbuf.begin(), recvbuf.size(), MPI_INT, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    if (_logger) _logger("MPI END: received chunk request from process " + std::to_string(recvbuf[0]));
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

    if (_logger)
        _logger("MPI BEGIN: serve chunk " + std::to_string(firstIndex) + ", " + std::to_string(numIndices)
                + " to process " + std::to_string(rank));
    std::array<size_t, 2> sendbuf{{firstIndex, numIndices}};
    MPI_Send(sendbuf.begin(), sendbuf.size(), MPI_UNSIGNED_LONG, rank, 1, MPI_COMM_WORLD);
    if (_logger) _logger("MPI END: served chunk");
#else
    (void)rank;
    (void)firstIndex;
    (void)numIndices;
    throwInvalidChunkInvocation();
#endif
}

//////////////////////////////////////////////////////////////////////

void ProcessManager::wait()
{
#ifdef BUILD_WITH_MPI
    if (isMultiProc())
    {
        if (_logger) _logger("MPI BEGIN: wait");
        MPI_Barrier(MPI_COMM_WORLD);
        if (_logger) _logger("MPI END: wait");
    }
#endif
}

//////////////////////////////////////////////////////////////////////

void ProcessManager::sumToAll(Array& arr)
{
#ifdef BUILD_WITH_MPI
    if (isMultiProc())
    {
        if (_logger) _logger("MPI BEGIN: sum to all of size " + std::to_string(arr.size()));
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
        if (_logger) _logger("MPI END: sum to all");
    }
#else
    (void)arr;
#endif
}

//////////////////////////////////////////////////////////////////////

void ProcessManager::sumToRoot(Array& arr, bool wait)
{
#ifdef BUILD_WITH_MPI
    if (isMultiProc())
    {
        if (_logger) _logger("MPI BEGIN: sum to root of size " + std::to_string(arr.size()));
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
        if (wait) MPI_Barrier(MPI_COMM_WORLD);
        if (_logger) _logger("MPI END: sum to root");
    }
#else
    (void)arr;
#endif
}

//////////////////////////////////////////////////////////////////////

void ProcessManager::broadcastAllToAll(std::function<void(vector<double>&)> producer,
                                       std::function<void(const vector<double>&)> consumer)
{
#ifdef BUILD_WITH_MPI
    if (isMultiProc())
    {
        if (_logger) _logger("MPI BEGIN: broadcast all to all");

        // allocate room for data to be sent and received
        vector<double> data;
        size_t datasize = 0;

        // iterate over all processes
        for (int k = 0; k != size(); ++k)
        {
            // if it is our turn to send, produce the data
            if (k == rank())
            {
                data.clear();
                producer(data);
                datasize = data.size();
            }

            // communicate the size of the data
            MPI_Bcast(&datasize, 1, MPI_UNSIGNED_LONG, k, MPI_COMM_WORLD);
            data.resize(datasize);

            // communicate the data itself, splitting it in maxMessageSize chunks if needed
            double* curdata = data.data();
            size_t remaining = datasize;
            while (remaining > maxMessageSize)
            {
                MPI_Bcast(curdata, maxMessageSize, MPI_DOUBLE, k, MPI_COMM_WORLD);
                remaining -= maxMessageSize;
                curdata += maxMessageSize;
            }
            if (remaining)
            {
                MPI_Bcast(curdata, remaining, MPI_DOUBLE, k, MPI_COMM_WORLD);
            }

            // unless it was our turn to send, consume the data
            if (k != rank())
            {
                consumer(data);
            }
        }
        if (_logger) _logger("MPI END: broadcast all to all");
    }
#else
    (void)producer;
    (void)consumer;
#endif
}

//////////////////////////////////////////////////////////////////////
