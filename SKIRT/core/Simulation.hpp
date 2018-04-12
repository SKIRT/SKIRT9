/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef SIMULATION_HPP
#define SIMULATION_HPP

#include "SimulationItem.hpp"
#include "ConsoleLog.hpp"
#include "FilePaths.hpp"
#include "ParallelFactory.hpp"
#include "PeerToPeerCommunicator.hpp"
#include "Random.hpp"
#include "Units.hpp"

////////////////////////////////////////////////////////////////////

/** Simulation is the abstract base class for a simulation item that represents a complete
    simulation and sits at the top of a run-time simulation hierarchy (i.e. it has no parent). A
    Simulation instance holds a number of essential simulation-wide property instances. Some of
    these (a random number generator and a system of units) are discoverable and hence fully
    user-configurable. The other properties (a file paths object, a logging mechanism, a parallel
    factory, and a peer-to-peer communicator) are not discoverable. When a Simulation instance is
    constructed, a default instance is created for each of these properties. A reference to these
    property instances can be retrieved through the corresponding getter, and in some cases, the
    property can be further configured under program control (e.g., to set the input and output
    file paths for the simulation).

    Specifically, when a Simulation instance is constructed, the \em filePaths property is set to
    an instance of the FilePaths class with default paths and no filename prefix; the \em log
    attribute is set to an instance of the Console class; the \em parallelFactory attribute is set
    to an instance of the ParallelFactory class with the default maximum number of parallel
    threads; and the \em communicator attribute is set to an instance of the PeerToPeerCommunicator
    class. */
class Simulation : public SimulationItem
{
    ITEM_ABSTRACT(Simulation, SimulationItem, "the simulation")

    PROPERTY_ITEM(random, Random, "the random number generator")
        ATTRIBUTE_DEFAULT_VALUE(random, "Random")

    PROPERTY_ITEM(units, Units, "the units system")
        ATTRIBUTE_DEFAULT_VALUE(units, "ExtragalacticUnits")

    ITEM_END()

    //======== Construction - Setup - Run - Destruction  ===========

public:
    /** This function performs setup for the complete simulation hierarchy. It invokes the setup()
        function defined in the SimulationItem base class, surrounded by start/finish log messages.
        It is recommended to use the setupAndRun() function rather than setup() and run() separately.
        */
    void setup();

    /** This function performs the simulation by invoking the runSelf() function to be defined in a
        subclass, surrounded by start/finish log messages. The setup() function must have been called
        before invoking run(). It is recommended to use the setupAndRun() function rather than
        setup() and run() separately. */
    void run();

    /** This function performs setup and executes the simulation by invoking setup() and run() in
        succession. */
    void setupAndRun();

protected:
    /** This function actually runs the simulation, assuming that setup() has been already performed.
        Its implementation must be provided by a subclass. */
    virtual void runSelf() = 0;

    //======== Getters for Non-Discoverable Attributes =======

public:
    /** Returns the PeerToPeerCommunicator of the simulation. */
    PeerToPeerCommunicator* communicator() const;

    /** Returns the logging mechanism for this simulation hierarchy. */
    Log* log() const;

    /** Returns the input/output file paths object for this simulation hierarchy. */
    FilePaths* filePaths() const;

    /** Returns the logging mechanism for this simulation hierarchy. */
    ParallelFactory* parallelFactory() const;

    //======================== Data Members ========================

private:
    // data members
    PeerToPeerCommunicator* _communicator{ new PeerToPeerCommunicator(this) };
    Log* _log{ new ConsoleLog(this) };
    FilePaths* _paths{ new FilePaths(this) };
    ParallelFactory* _factory{ new ParallelFactory(this) };
};

////////////////////////////////////////////////////////////////////

#endif
