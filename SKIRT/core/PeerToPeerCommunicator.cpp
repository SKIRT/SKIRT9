/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "PeerToPeerCommunicator.hpp"
#include "Log.hpp"
#include "ProcessManager.hpp"

////////////////////////////////////////////////////////////////////

namespace
{
    const int ROOT = 0;
}

////////////////////////////////////////////////////////////////////

PeerToPeerCommunicator::PeerToPeerCommunicator(SimulationItem* parent)
{
    parent->addChild(this);
}

////////////////////////////////////////////////////////////////////

void PeerToPeerCommunicator::setDataParallel(bool dataParallel)
{
    _dataParallel = dataParallel;
}

////////////////////////////////////////////////////////////////////

void PeerToPeerCommunicator::sum(Array& arr)
{
    if (!isMultiProc()) return;

    ProcessManager::sum(&(arr[0]),arr.size(),0);
}

////////////////////////////////////////////////////////////////////

void PeerToPeerCommunicator::sumAll(Array& arr)
{
    if (!isMultiProc()) return;

    ProcessManager::sumAll(&(arr[0]),arr.size());
}

////////////////////////////////////////////////////////////////////

void PeerToPeerCommunicator::sumAll(double& dbl)
{
    if (!isMultiProc()) return;

    ProcessManager::sumAll(&dbl,1);
}

void PeerToPeerCommunicator::orAll(bool& b)
{
    if (!isMultiProc()) return;

    ProcessManager::orAll(&b);
}

////////////////////////////////////////////////////////////////////

void PeerToPeerCommunicator::broadcast(Array& arr, int sender)
{
    if (!isMultiProc()) return;

    ProcessManager::broadcast(&(arr[0]),arr.size(),sender);
}

////////////////////////////////////////////////////////////////////

void PeerToPeerCommunicator::broadcast(int& value, int sender)
{
    if (!isMultiProc()) return;

    ProcessManager::broadcast(&value,sender);
}

////////////////////////////////////////////////////////////////////

void PeerToPeerCommunicator::gatherWithPattern(const double* sendBuffer, size_t sendCount,
                                     double* recvBuffer, int recvRank, size_t recvLength,
                                     const std::vector<std::vector<int>>& recvDisplacements)
{
    if (!isMultiProc()) return;

    ProcessManager::gatherWithPattern(sendBuffer, sendCount, recvBuffer, recvRank, recvLength, recvDisplacements);
}

////////////////////////////////////////////////////////////////////

void PeerToPeerCommunicator::displacedBlocksAllToAll(const double* sendBuffer, size_t sendCount, size_t sendLength,
                                                     std::vector<std::vector<int>>& sendDisplacements,
                                                     size_t sendExtent,
                                                     double* recvBuffer, size_t recvCount, size_t recvLength,
                                                     std::vector<std::vector<int>>& recvDisplacements,
                                                     size_t recvExtent)
{
    if (!isMultiProc()) return;

    ProcessManager::displacedBlocksAllToAll(sendBuffer, sendCount, sendLength, sendDisplacements, sendExtent,
                                            recvBuffer, recvCount, recvLength, recvDisplacements, recvExtent);
}

////////////////////////////////////////////////////////////////////

int PeerToPeerCommunicator::root()
{
    return ROOT;
}

////////////////////////////////////////////////////////////////////

bool PeerToPeerCommunicator::isRoot()
{
    return (rank() == ROOT);
}

////////////////////////////////////////////////////////////////////

void PeerToPeerCommunicator::wait(string scope)
{
    if (!isMultiProc()) return;

    find<Log>()->info("Waiting for other processes to finish " + scope + "...");

    ProcessManager::barrier();
}

////////////////////////////////////////////////////////////////////

bool PeerToPeerCommunicator::dataParallel()
{
    return _dataParallel;
}

////////////////////////////////////////////////////////////////////
