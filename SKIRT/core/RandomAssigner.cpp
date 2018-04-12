/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "RandomAssigner.hpp"
#include "FatalError.hpp"
#include "PeerToPeerCommunicator.hpp"
#include "Random.hpp"
#include "SequentialAssigner.hpp"

////////////////////////////////////////////////////////////////////

RandomAssigner::RandomAssigner(size_t size, SimulationItem* parent)
    : ProcessAssigner(size, parent)
{
    if (!communicator()) throw FATALERROR("Could not find an object of type PeerToPeerCommunicator in the simulation hierarchy");
    _random = find<Random>();

    _assignment.clear();
    _assignment.resize(size);
    _values.clear();
    _values.reserve(static_cast<size_t>(1.2*size/communicator()->size()));

    // For each value in a certain subset of 'size', let this process determine a random process rank
    SequentialAssigner* helpassigner = new SequentialAssigner(size, this);
    for (size_t i = 0; i < helpassigner->assigned(); i++)
    {
        int rank = static_cast<int>(floor(_random->uniform() * communicator()->size()));
        _assignment[helpassigner->absoluteIndex(i)] = rank;
    }

    // Communication of the randomly determined ranks
    for (size_t j = 0; j < size; j++)
    {
        int sender = helpassigner->rankForIndex(j);
        communicator()->broadcast(_assignment[j], sender);

        // If the process assigned to this value is this process, add the value to the list
        if (_assignment[j] == communicator()->rank())
        {
            _values.push_back(j);
        }
    }

    // Set the number of values assigned to this process
    setAssigned(_values.size());
}

////////////////////////////////////////////////////////////////////

size_t RandomAssigner::absoluteIndex(size_t relativeIndex) const
{
    return _values[relativeIndex];
}

////////////////////////////////////////////////////////////////////

size_t RandomAssigner::relativeIndex(size_t absoluteIndex) const
{
    return std::find(_values.begin(), _values.end(), absoluteIndex) - _values.begin();
}

////////////////////////////////////////////////////////////////////

int RandomAssigner::rankForIndex(size_t index) const
{
    return _assignment[index];
}

////////////////////////////////////////////////////////////////////
