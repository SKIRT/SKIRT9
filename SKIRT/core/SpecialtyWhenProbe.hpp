/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef SPECIALTYWHENPROBE_HPP
#define SPECIALTYWHENPROBE_HPP

#include "SpecialtyProbe.hpp"

////////////////////////////////////////////////////////////////////

/** SpecialtyWhenProbe is a base class for specialty probes that allow the user to decide whether
    the probe should be performed after setup or after the full simulation run. */
class SpecialtyWhenProbe : public SpecialtyProbe
{
    /** The enumeration type indicating when probing occurs. */
    ENUM_DEF(ProbeAfter, Setup, Run)
        ENUM_VAL(ProbeAfter, Setup, "after setup")
        ENUM_VAL(ProbeAfter, Run, "after the complete simulation run")
    ENUM_END()

    ITEM_ABSTRACT(SpecialtyWhenProbe, SpecialtyProbe, "a specialty when probe")

        ATTRIBUTE_SUB_PROPERTIES_HERE(SpecialtyWhenProbe)

        PROPERTY_ENUM(probeAfter, ProbeAfter, "perform the probe after")
        ATTRIBUTE_DEFAULT_VALUE(probeAfter, "Setup")
        ATTRIBUTE_DISPLAYED_IF(probeAfter, "DynamicState")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function returns an enumeration indicating when probing for this probe should be
        performed corresponding to the configured value of the \em probeAfter property. */
    When when() const override;
};

////////////////////////////////////////////////////////////////////

#endif
