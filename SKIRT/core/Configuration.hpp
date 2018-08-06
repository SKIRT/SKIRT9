/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef CONFIGURATION_HPP
#define CONFIGURATION_HPP

#include "SimulationItem.hpp"

////////////////////////////////////////////////////////////////////

/** Configuration is a helper class that serves as a central clearing house for overall simulation
    configuration options, including those offered by all members of the SimulationMode class
    hierarchy.

    Each MonteCarloSimulation holds a single Configuration object. During setup, it retrieves many
    properties and options from the simulation hierarchy, verifying consistency of the
    configuration and flagging any conflicts while doing so. Once this process has completed, the
    Configuration object offers getters for these retrieved properties to any of the other
    simulation items in the hierarchy. The setup() function of the Configuration object is invoked
    at the very early stages of simulation setup, so that it is safe for other simulation items to
    retrieve information from the Configuration object during setup.

    The Configuration class is based on SimulationItem so that it can be part of a simulation item
    hierarchy, however it is not discoverable because it is not intended to be selected or
    configured by the user. */
class Configuration : public SimulationItem
{
    //============= Construction - Setup - Destruction =============

public:
    /** This constructor creates a Configuration object that is hooked up as a child to the
        specified parent in the simulation hierarchy, so that it will automatically be deleted. The
        setup() function is \em not called by this constructor. */
    explicit Configuration(SimulationItem* parent);

protected:
    /** This function retrieves properties and options from the simulation hierarchy and stores the
        resulting values internally so that they can be returned by any of the getters with minimal
        overhead. During this process, the function also verifies the consistency of the simulation
        configuration, for example checking the configuration against the requirements of the
        selected SimulationMode subclass. If any conflicts are found, the function throws a fatal
        error. */
    void setupSelfBefore() override;

    //======== Setters that override the user configuration =======

public:
    /** This function puts the simulation in emulation mode. Specifically, it sets a flag that can
        be queried by other simulation items, it sets the number of photon packets to zero, and if
        iteration over the simulation state is enabled, it forces the number of iterations to one.
        */
    void setEmulationMode();

    //=========== Getters for configuration properties ============

public:
    /** Returns true if the simulation has been put in emulation mode. */
    bool emulationMode() { return _emulationMode; }

    /** Returns the number of photon packets launched per primary emission simulation segment. */
    double numPrimaryPackets() const { return _numPrimaryPackets; }

    /** Returns the minimum weight reduction factor before a photon packet is terminated. */
    double minWeightReduction() const { return _minWeightReduction; }

    /** Returns the minimum number of forced scattering events before a photon packet is
        terminated. */
    int minScattEvents() const { return _minScattEvents; }

    /** Returns the fraction of path lengths sampled from a linear rather than an exponential
        distribution. */
    double pathLengthBias() const { return _pathLengthBias; }

    /** Returns true if the media in the simulation support polarization. */
    bool hasPolarization() const { return _hasPolarization; }

    //======================== Data Members ========================

private:
    // general
    bool _emulationMode{false};

    // no media
    double _numPrimaryPackets{0.};

    // extinction only
    double _minWeightReduction{1e4};
    int _minScattEvents{0};
    double _pathLengthBias{0.5};
    bool _hasPolarization{false};
};

////////////////////////////////////////////////////////////////////

#endif
