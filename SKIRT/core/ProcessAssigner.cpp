/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "ProcessAssigner.hpp"
#include "PeerToPeerCommunicator.hpp"

////////////////////////////////////////////////////////////////////

ProcessAssigner::ProcessAssigner(size_t size, SimulationItem* parent)
    : _total(size)
{
    parent->addChild(this);
    setup();
}

////////////////////////////////////////////////////////////////////

void ProcessAssigner::setupSelfBefore()
{
    SimulationItem::setupSelfBefore();

    // cache a pointer to the PeerToPeerCommunicator
    _comm = find<PeerToPeerCommunicator>();
}

////////////////////////////////////////////////////////////////////

void ProcessAssigner::setAssigned(size_t assigned)
{
    _assigned = assigned;
}

////////////////////////////////////////////////////////////////////

size_t ProcessAssigner::total() const
{
    return _total;
}

////////////////////////////////////////////////////////////////////

size_t ProcessAssigner::assigned() const
{
    return _assigned;
}

////////////////////////////////////////////////////////////////////

size_t ProcessAssigner::assignedForRank(int rank) const
{
    size_t result = 0;

    for (size_t absoluteIndex=0; absoluteIndex<total(); absoluteIndex++)
        if (rank == rankForIndex(absoluteIndex)) result++;
    return result;
}

////////////////////////////////////////////////////////////////////

bool ProcessAssigner::validIndex(size_t absoluteIndex) const
{
    return _comm->rank() == rankForIndex(absoluteIndex);
}

////////////////////////////////////////////////////////////////////

vector<int> ProcessAssigner::indicesForRank(int rank) const
{
    vector<int> result;
    result.reserve(_assigned);

    // add all the absolute indices that correspond to the given rank
    for (size_t absoluteIndex=0; absoluteIndex<total(); absoluteIndex++)
        if (rank == rankForIndex(absoluteIndex))
            result.push_back(absoluteIndex);

    return result;
}

////////////////////////////////////////////////////////////////////
