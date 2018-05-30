/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef MONTECARLOSIMULATION_HPP
#define MONTECARLOSIMULATION_HPP

#include "Simulation.hpp"
#include "InstrumentSystem.hpp"
#include "ProbeSystem.hpp"
#include "SourceSystem.hpp"
#include <atomic>

//////////////////////////////////////////////////////////////////////

/** The MonteCarloSimulation class is the general abstract base class describing Monte Carlo
    simulations. Running a Monte Carlo simulation with SKIRT essentially comes down to constructing
    an instance of one of the subclasses of the MonteCarloSimulation base class and invoking the
    setupAndRun() function on it. The MonteCarloSimulation class ... . */
class MonteCarloSimulation : public Simulation
{
    ITEM_CONCRETE(MonteCarloSimulation, Simulation, "a Monte Carlo simulation")

    ATTRIBUTE_SUB_PROPERTIES_HERE(MonteCarloSimulation)

    PROPERTY_ITEM(sourceSystem, SourceSystem, "the source system")
        ATTRIBUTE_DEFAULT_VALUE(sourceSystem, "SourceSystem")

    PROPERTY_ITEM(instrumentSystem, InstrumentSystem, "the instrument system")
        ATTRIBUTE_DEFAULT_VALUE(instrumentSystem, "InstrumentSystem")

    PROPERTY_ITEM(probeSystem, ProbeSystem, "the probe system")
        ATTRIBUTE_DEFAULT_VALUE(probeSystem, "ProbeSystem")

    PROPERTY_DOUBLE(numPackets, "the number of photon packets launched per simulation segment")
        ATTRIBUTE_MIN_VALUE(numPackets, "[0")
        ATTRIBUTE_MAX_VALUE(numPackets, "1e19]")
        ATTRIBUTE_DEFAULT_VALUE(numPackets, "1e6")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function performs setup for the complete simulation hierarchy. It calls the regular
        setup() function and notifies the probe system when setup has been completed. */
    void setupSimulation() override;

    /** This function performs setup for the MonteCarloSimulation object; it caches some frequently
        used pointers. */
    void setupSelfBefore() override;

    //======== Setters & Getters for Discoverable Attributes =======

    /** \fn numPackets
        The number of photon packets is specified as a double-precision floating point number
        rather than as a 64-bit integer to avoid implementing yet another discoverable property
        type. As a side benefit, one can use exponential notation to specify a large number of
        photon packets. Also, note that a double can exactly represent all integers up to 9e15. The
        maximum number of photon packets is somewhat arbitrarily set to 1e19 because that number is
        close to the maximum number representable with a 64-bit unsigned integer. */

public:
    /** This function puts the simulation in emulation mode. Specifically, it sets an internal flag
        that can be queried other classes and it sets the number of photon packets to zero. */
    void setEmulationMode();

    /** This function returns true if the simulation has been put in emulation mode. */
    bool emulationMode();

    //======================== Other Functions =======================

public:
    /** This function returns the dimension of the simulation, which depends on the (lack of)
        symmetry in the media geometries. A value of 1 means spherical symmetry, 2 means axial
        symmetry and 3 means none of these symmetries. The media component with the least symmetry
        (i.e. the highest dimension) determines the result for the whole simulation. */
    int dimension() const;

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
    void initProgress(string segment, size_t numPackets);

    /** This function logs a progress message for the segment specified in the initprogress()
        function if the previous message was issued at least 3 seconds ago. The function
        must be called regularly while processing photon packets. The argument specifies the
        number of photon packets processed since the most recent invocation in the same thread. */
    void logProgress(size_t extraDone);

    /** This function launches the specified chunk of photon packets from primary sources. */
    void doPrimaryEmissionChunk(size_t firstIndex, size_t numIndices);

    //======================== Data Members ========================

private:
    // *** data member to remember whether emulation mode is enabled
    bool _emulationMode{false};

    // *** data members used by the XXXprogress() functions in this class ***
    string _segment;               // a string identifying the photon shooting segment for use in the log message
    size_t _numTotal;              // the total number of photon packages to be processed for this segment
    std::atomic<size_t> _numDone;  // the number of photon packages processed so far (for all threads)
};

////////////////////////////////////////////////////////////////////

#endif
