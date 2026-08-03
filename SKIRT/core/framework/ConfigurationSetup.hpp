/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef CONFIGURATIONSETUP_HPP
#define CONFIGURATIONSETUP_HPP

#include "Configuration.hpp"

////////////////////////////////////////////////////////////////////

/** ConfigurationSetup is the concrete subclass of Configuration that performs the actual
    retrieval of overall simulation configuration options.

    Each MonteCarloSimulation holds a single ConfigurationSetup object, accessed through the
    Configuration interface. During setup, it retrieves many properties and options from the
    simulation hierarchy, verifying consistency of the configuration and flagging any conflicts
    while doing so, and stores the results in the (protected) data members inherited from
    Configuration so that they can be returned by any of the getters declared there with minimal
    overhead. The setup() function of this object is invoked at the very early stages of
    simulation setup, so that it is safe for other simulation items to retrieve information from
    the Configuration object during their own setup. */
class ConfigurationSetup : public Configuration
{
    //============= Construction - Setup - Destruction =============

public:
    /** This constructor creates a ConfigurationSetup object that is hooked up as a child to the
        specified parent in the simulation hierarchy, so that it will automatically be deleted. The
        setup() function is \em not called by this constructor. */
    explicit ConfigurationSetup(SimulationItem* parent);

protected:
    /** This function retrieves properties and options from the simulation hierarchy and stores the
        resulting values internally so that they can be returned by any of the getters with minimal
        overhead. During this process, the function also verifies the consistency of the simulation
        configuration, for example checking the configuration against the requirements of the
        selected simulation mode. If any conflicts are found, the function throws a fatal error. */
    void setupSelfBefore() override;

    /** This function logs some aspects of the configuration as information to the user. */
    void setupSelfAfter() override;

    //=========== Getters requiring access to the live hierarchy ===========

public:
    /** Returns a wavelength range that covers all wavelengths possibly used in the simulation for
        photon transport or for otherwise probing material properties (e.g. optical depth). See the
        base class for the full description. */
    Range simulationWavelengthRange() const override;

    /** Returns a list of wavelengths that are explicitly or indirectly mentioned by the simulation
        configuration. See the base class for the full description. */
    vector<double> simulationWavelengths() const override;
};

////////////////////////////////////////////////////////////////////

#endif
