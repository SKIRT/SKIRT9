/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef MONTECARLOSIMULATION_HPP
#define MONTECARLOSIMULATION_HPP

#include "Simulation.hpp"
#include "Configuration.hpp"
#include "InstrumentSystem.hpp"
#include "MediumSystem.hpp"
#include "ProbeSystem.hpp"
#include "SimulationMode.hpp"
#include "SourceSystem.hpp"
#include <atomic>

//////////////////////////////////////////////////////////////////////

/** The MonteCarloSimulation class is the top-level class describing a SKIRT simulation. Running a
    Monte Carlo simulation with SKIRT essentially comes down to constructing an instance of the
    MonteCarloSimulation class and invoking the setupAndRun() function on it.

    The MonteCarloSimulation class holds the source, media, instrument and probe systems,
    implements the core aspects of the photon packet life-cycle, and manages the iterative
    processes in the simulation (including phases, iterations and segments).

    TODO: more documentation ...

    The MonteCarloSimulation class also holds the non-discoverable \em config property, which is
    automatically set to an instance of the Configuration class. The setup() function of the config
    object is invoked at the very early stages of overall simulation setup, so that it can
    initialize its internal state to reflect the simulation configuration. As a result, it is safe
    for other simulation items to retrieve information from the config object during setup. */
class MonteCarloSimulation : public Simulation
{
    ITEM_CONCRETE(MonteCarloSimulation, Simulation, "a Monte Carlo simulation")

    PROPERTY_ITEM(mode, SimulationMode, "the overall simulation mode")
        ATTRIBUTE_DEFAULT_VALUE(mode, "ExtinctionOnlyMode")

    PROPERTY_ITEM(sourceSystem, SourceSystem, "the source system")
        ATTRIBUTE_DEFAULT_VALUE(sourceSystem, "SourceSystem")

    PROPERTY_ITEM(mediumSystem, MediumSystem, "the medium system")
        ATTRIBUTE_DEFAULT_VALUE(mediumSystem, "MediumSystem")
        ATTRIBUTE_OPTIONAL(mediumSystem)

    PROPERTY_ITEM(instrumentSystem, InstrumentSystem, "the instrument system")
        ATTRIBUTE_DEFAULT_VALUE(instrumentSystem, "InstrumentSystem")

    PROPERTY_ITEM(probeSystem, ProbeSystem, "the probe system")
        ATTRIBUTE_DEFAULT_VALUE(probeSystem, "ProbeSystem")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function performs setup for the complete simulation hierarchy. It calls the regular
        setup() function and notifies the probe system when setup has been completed. */
    void setupSimulation() override;

    /** This function performs initial setup for the MonteCarloSimulation object; it caches some
        frequently used pointers. */
    void setupSelfBefore() override;

    /** This function performs final setup for the MonteCarloSimulation object; it logs the
        dimension of the simulation. */
    void setupSelfAfter() override;

    //======== Getters for Non-Discoverable Properties =======

public:
    /** Returns the Configuration object for this simulation hierarchy. */
    Configuration* config() const;

    //======================== Other Functions =======================

protected:
    /** This function actually runs the simulation, assuming setup has been completed. */
    void runSimulation() override;

private:
    /** In a multi-processing environment, this function logs a message and waits for all processes
        to finish the work (i.e. it places a barrier). The string argument is included in the log
        message to indicate the scope of work that is being finished. If there is only a single
        process, the function does nothing. */
    void wait(string scope);

    /** This function initializes the progress counter used in logprogress() for the specified
        segment and logs the number of photon packets to be processed. */
    void initProgress(string segment, size_t numTotal);

    /** This function logs a progress message for the segment specified in the initprogress()
        function if the previous message was issued at least 3 seconds ago. The function
        must be called regularly while processing photon packets. The argument specifies the
        number of photon packets processed so far. */
    void logProgress(size_t numDone);

    /** This function launches the specified chunk of photon packets from primary sources. */
    void doPrimaryEmissionChunk(size_t firstIndex, size_t numIndices);

    //======================== Data Members ========================

private:
    // non-discoverable simulation items
    Configuration* _config{ new Configuration(this) };

    // data members used by the XXXprogress() functions in this class
    string _segment;               // a string identifying the photon shooting segment for use in the log message
    size_t _numTotal;              // the total number of photon packages to be processed for this segment
};

////////////////////////////////////////////////////////////////////

#endif
