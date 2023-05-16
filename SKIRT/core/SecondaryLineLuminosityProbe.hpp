/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef SECONDARYLINELUMINOSITYPROBE_HPP
#define SECONDARYLINELUMINOSITYPROBE_HPP

#include "SpatialGridFormProbe.hpp"

////////////////////////////////////////////////////////////////////

/** SecondaryLineLuminosityProbe probes the line luminosities of secondary sources as discretized
    on the spatial grid of the simulation. It produces a separate output file for each medium
    component with an associated material mix that supports secondary line emission.

    For each relevant medium component, the probe outputs a compound quantity representing the
    intrinsic bolometric, line-profile-integrated luminosities for each of the lines from which the
    medium emits, ignoring any kinematic effects. These luminosities are listed with the
    corresponding central line wavelengths, frequencies, or energies, depending on the \em
    wavelengthOutputStyle configured for the simulation.

    The probe can be used with any Form subclass. When associated with a form that samples the
    luminosity in a spatial grid cell or at a set of positions, such as for a linear or planar cut,
    the probe outputs a bolometric luminosity volume density (with SI units of W/m3). When
    associated with a form that projects the luminosity along a path, the probe outputs a
    bolometric luminosity surface (or column) density (with SI units of W/m2). */
class SecondaryLineLuminosityProbe : public SpatialGridFormProbe
{
    /** The enumeration type indicating when probing occurs. */
    ENUM_DEF(ProbeAfter, Run, Secondary)
        ENUM_VAL(ProbeAfter, Run, "after the complete simulation run")
        ENUM_VAL(ProbeAfter, Secondary, "after each iteration over secondary emission")
    ENUM_END()

    ITEM_CONCRETE(SecondaryLineLuminosityProbe, SpatialGridFormProbe,
                  "internal spatial grid: secondary line luminosity")
        ATTRIBUTE_TYPE_DISPLAYED_IF(SecondaryLineLuminosityProbe, "Level2&SpatialGrid&GasEmission")

        PROPERTY_ENUM(probeAfter, ProbeAfter, "perform the probe after")
        ATTRIBUTE_DEFAULT_VALUE(probeAfter, "Run")
        ATTRIBUTE_DISPLAYED_IF(probeAfter, "IterateSecondary")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function returns an enumeration indicating when probing for this probe should be
        performed corresponding to the configured value of the \em probeAfter property. */
    When when() const override;

    //======================== Other Functions =======================

public:
    /** This function performs probing. */
    void probe() override;
};

////////////////////////////////////////////////////////////////////

#endif
