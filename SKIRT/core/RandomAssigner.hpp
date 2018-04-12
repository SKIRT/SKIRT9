/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef RANDOMASSIGNER_HPP
#define RANDOMASSIGNER_HPP

#include "ProcessAssigner.hpp"
class Random;

//////////////////////////////////////////////////////////////////////

/** The RandomAssigner class is a subclass of the ProcessAssigner class, representing objects that
    assign work to different processes. An object of the RandomAssigner class distributes the work
    amongst the different processes in a random way. For each value (part of work), the assigned
    process is determined by drawing a uniform random number. Because the state of the random
    number generator can be different in each process, each process performs the random assignment
    for a subset of the different values, after which the results are gathered across all
    processes. */
class RandomAssigner : public ProcessAssigner
{
    //============= Construction - Setup - Destruction =============

public:
    /** This constructor can be invoked by SKIRT classes that wish to hard-code the creation of a
        new ProcessAssigner object of this type. Before the constructor returns, the newly created
        object is hooked up as a child to the specified parent in the simulation hierarchy (so it
        will automatically be deleted) and the setup of the ProcessAssigner base class is invoked,
        which sets the _comm attribute that points to the object of type PeerToPeerCommunicator
        that is found in the simulation hierarchy. The second argument takes the number of parts of
        work \f$n\f$ that need to be performed. For each part of work \f$i\f$ (with \f$i\f$ ranging
        from 0 to \f$n-1\f$), the rank of the assigned process is determined, by drawing random
        numbers from a uniform distribution. By using a SequentialAssigner, the range of indices
        \f$i\f$ is split in different domains, so that each process must draw these random numbers
        only for one particular domain. Subsequently, the randomly determined ranks are
        synchronized across all processes, so that every process knows which process is assigned to
        any part of the work (stored in the _assignment vector). Each process then determines which
        values it is assigned to itself, and stores these values in a seperate list (the _values
        attribute). */
    explicit RandomAssigner(size_t size, SimulationItem* parent);

    //======================== Other Functions =======================

public:
    /** This function takes the relative index of a certain part of the work assigned to this process
        as an argument and returns the absolute index of that part, a value from zero to the total
        amount of parts that need to be executed in the simulation. This is done by simply looking up
        the value in the _values list. */
    size_t absoluteIndex(size_t relativeIndex) const override;

    /** This function takes the absolute index of a certain part of the work as an argument and returns
        the relative index of that part, a value from zero to the number of parts that were assigned to
        this process. This is done by performing a search through the _values list for the value that
        matches the absoluteIndex. The relative index that is returned corresponds to the index of that
        value in the _values list. */
    size_t relativeIndex(size_t absoluteIndex) const override;

    /** This function returns the rank of the process that is assigned to a certain part of the work.
        This part is identified by its absolute index, passed to this function as an argument. This is
        done by simply looking up the corresponding rank in the _assignment vector. */
    int rankForIndex(size_t index) const override;

    //======================== Data Members ========================

private:
    Random* _random{nullptr};  // a pointer to the Random object of the simulation
    vector<int> _assignment;   // for each value, this vector defines the rank of the assigned process
    vector<size_t> _values;    // a list of the values assigned to this process
};

//////////////////////////////////////////////////////////////////////

#endif
