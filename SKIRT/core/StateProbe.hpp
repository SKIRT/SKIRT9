/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef STATEPROBE_HPP
#define STATEPROBE_HPP

#include "Probe.hpp"

////////////////////////////////////////////////////////////////////

/** StateProbe is an abstract Probe subclass representing probes that query the medium state and
    thus, if the simulation has a dynamic medium state, may produce a different result depending on
    when the probe is performed. Accordingly, this class offers an option for the user to decide
    when the probe should be performed (after setup or after the full simulation run), and invokes
    the subclass probe() function at the user-configured moment. */
class StateProbe : public Probe
{
    /** The enumeration type indicating when probing occurs. */
    ENUM_DEF(ProbeAfter, Setup, Run)
        ENUM_VAL(ProbeAfter, Setup, "after setup")
        ENUM_VAL(ProbeAfter, Run, "after the complete simulation run")
    ENUM_END()

    ITEM_ABSTRACT(StateProbe, Probe, "a probe that queries the medium state")

        ATTRIBUTE_SUB_PROPERTIES_HERE()

        PROPERTY_ENUM(probeAfter, ProbeAfter, "perform the probe after")
        ATTRIBUTE_DEFAULT_VALUE(probeAfter, "Setup")
        ATTRIBUTE_DISPLAYED_IF(probeAfter, "DynamicState")

    ITEM_END()

    //======================== Other Functions =======================

public:
    /** This function performs probing after setup. It invokes probe() only if the \em probeAfter
        property is set to Setup. */
    void probeSetup() override final;

    /** This function performs probing after all photon packets have been emitted and detected. It
        invokes probe() only if the \em probeAfter property is set to Run. */
    void probeRun() override final;

protected:
    /** This function must be implemented in subclasses to perform the probing; it is called from
        probeSetup() or probeRun() depending on the value of the \em probeAfter property. */
    virtual void probe() = 0;
};

////////////////////////////////////////////////////////////////////

#endif
