/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef SIMULATION_HPP
#define SIMULATION_HPP

#include "ConsoleLog.hpp"
#include "FilePaths.hpp"
#include "ParallelFactory.hpp"
#include "Random.hpp"
#include "SimulationItem.hpp"
#include "Units.hpp"

////////////////////////////////////////////////////////////////////

/** Simulation is the abstract base class for a simulation item that represents a complete
    simulation and sits at the top of a run-time simulation hierarchy (i.e. it has no parent). A
    Simulation instance holds a number of essential simulation-wide property instances. Some of
    these (a random number generator and a system of units) are discoverable and hence fully
    user-configurable. The other properties (a file paths object, a logging mechanism, and a
    parallel factory) are not discoverable. When a Simulation instance is constructed, a default
    instance is created for each of these properties. A reference to these property instances can
    be retrieved through the corresponding getter, and in some cases, the property can be further
    configured under program control (e.g., to set the input and output file paths for the
    simulation).

    Specifically, when a Simulation instance is constructed, the \em log property is set to an
    instance of the ConsoleLog class; the \em filePaths property is set to an instance of the
    FilePaths class with default paths and no filename prefix; and the \em parallelFactory property
    is set to an instance of the ParallelFactory class with the default maximum number of parallel
    threads. */
class Simulation : public SimulationItem
{
    /** The enumeration type indicating the user experience level:
        - Level 1: Basic
    */
    ENUM_DEF(UserLevel, Basic, Regular, Expert)
        ENUM_VAL(UserLevel, Basic, "Basic: for beginning users (hides many options)")
        ENUM_VAL(UserLevel, Regular, "Regular: for regular users (hides esoteric options)")
        ENUM_VAL(UserLevel, Expert, "Expert: for expert users (hides no options)")
    ENUM_END()

    ITEM_ABSTRACT(Simulation, SimulationItem, "the simulation")

        PROPERTY_ENUM(userLevel, UserLevel, "the user experience level")
        ATTRIBUTE_DEFAULT_VALUE(userLevel, "Regular")
        ATTRIBUTE_INSERT(userLevel,
                         "userLevelBasic:Level1;userLevelRegular:Level1,Level2;userLevelExpert:Level1,Level2,Level3")

        PROPERTY_ITEM(random, Random, "the random number generator")
        ATTRIBUTE_DEFAULT_VALUE(random, "Random")
        ATTRIBUTE_DISPLAYED_IF(random, "Level3")

        PROPERTY_ITEM(units, Units, "the units system")
        ATTRIBUTE_DEFAULT_VALUE(units, "ExtragalacticUnits")

    ITEM_END()

    //======== Construction - Setup - Run - Destruction  ===========

public:
    /** This function performs setup for the complete simulation hierarchy by invoking the setup()
        function defined in the SimulationItem base class, and then runs the simulation by invoking
        the run() function which must be defined in a subclass. The complete operation is
        surrounded by start/finish log messages.

        It is highly recommended for the creator/manager of a simulation hierarchy to immediately
        call setupAndRun() on the Simulation instance rather than first calling the setup()
        function. */
    void setupAndRun();

protected:
    /** This function actually performs setup for the complete simulation hierarchy. Its
        implementation must be provided by a subclass. */
    virtual void setupSimulation() = 0;

    /** This function actually runs the simulation, assuming that setup has been already performed.
        Its implementation must be provided by a subclass. */
    virtual void runSimulation() = 0;

    //======== Getters for Non-Discoverable Attributes =======

public:
    /** Returns the logging mechanism for this simulation hierarchy. */
    Log* log() const;

    /** Returns the input/output file paths object for this simulation hierarchy. */
    FilePaths* filePaths() const;

    /** Returns the logging mechanism for this simulation hierarchy. */
    ParallelFactory* parallelFactory() const;

    //======================== Data Members ========================

private:
    // data members
    Log* _log{new ConsoleLog(this)};
    FilePaths* _paths{new FilePaths(this)};
    ParallelFactory* _factory{new ParallelFactory(this)};
};

////////////////////////////////////////////////////////////////////

#endif
