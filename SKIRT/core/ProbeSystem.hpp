/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef PROBESYSTEM_HPP
#define PROBESYSTEM_HPP

#include "Probe.hpp"

//////////////////////////////////////////////////////////////////////

/** A ProbeSystem instance maintains a list of zero or more probes and helps ensure that these probes are
invoked at the appropriate times. */
class ProbeSystem : public SimulationItem
{
    ITEM_CONCRETE(ProbeSystem, SimulationItem, "a probe system")

    PROPERTY_ITEM_LIST(probes, Probe, "the probes")
        ATTRIBUTE_OPTIONAL(instruments)

    ITEM_END()

    //======================== Other Functions =======================

public:
    /** This function invokes the function of the same name on all probes in the probe system. */
    void probeSetup();

    /** This function invokes the function of the same name on all probes in the probe system. */
    void probeSimulation();
};

////////////////////////////////////////////////////////////////////

#endif
