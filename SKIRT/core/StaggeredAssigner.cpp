/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "StaggeredAssigner.hpp"
#include "FatalError.hpp"
#include "PeerToPeerCommunicator.hpp"

////////////////////////////////////////////////////////////////////

StaggeredAssigner::StaggeredAssigner(size_t size, SimulationItem* parent)
    : ProcessAssigner(size, parent)
{
    if (!communicator()) throw FATALERROR("Could not find an object of type PeerToPeerCommunicator in the simulation hierarchy");

    size_t assigned = 0;
    for (size_t i = 0; i < size; i++)
    {
        int rank = i % communicator()->size();
        if (rank == communicator()->rank()) assigned++;
    }

    // Set the total number of assigned values
    setAssigned(assigned);
}

////////////////////////////////////////////////////////////////////

size_t StaggeredAssigner::absoluteIndex(size_t relativeIndex) const
{
    return communicator()->rank() + relativeIndex*communicator()->size();
}

////////////////////////////////////////////////////////////////////

size_t StaggeredAssigner::relativeIndex(size_t absoluteIndex) const
{
    return ((absoluteIndex - communicator()->rank()) / communicator()->size());
}

////////////////////////////////////////////////////////////////////

int StaggeredAssigner::rankForIndex(size_t index) const
{
    return (index % communicator()->size());
}

////////////////////////////////////////////////////////////////////
