/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef PROBE_HPP
#define PROBE_HPP

#include "SimulationItem.hpp"

////////////////////////////////////////////////////////////////////

/** Probe is an abstract class representing probes that output stuff. */
class Probe : public SimulationItem
{
    ITEM_ABSTRACT(Probe, SimulationItem, "a probe")
    PROPERTY_STRING(probeName, "the name for this probe")
    ITEM_END()

    //======================== Other Functions =======================

public:
    /** This function returns the probe name as human-readable name for the simulation item, so
        that it can be used in log messages to identify the probe and differentiate it from other
        probes. */
    string itemName() const override;

    /** This function performs probing after setup. */
    virtual void probeSetup();

    /** This function performs probing after the simulation. */
    virtual void probeRun();
};

////////////////////////////////////////////////////////////////////

#endif
