/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef SECONDARYDUSTLUMINOSITYPROBE_HPP
#define SECONDARYDUSTLUMINOSITYPROBE_HPP

#include "SpatialGridFormProbe.hpp"

////////////////////////////////////////////////////////////////////

/** SecondaryDustLuminosityProbe probes the specific luminosity aggregated over all dust medium
    components in the simulation at the wavelengths specified by the dust emission wavelength grid
    configured for the simulation. It produces output only if the simulation has at least one dust
    medium component and includes dust emission.

    When associated with a form that samples the luminosity per spatial cell or at a set of
    positions, such as for a linear or planar cut, the probe outputs a monochromatic luminosity
    volume density (with SI units of W/m/m3 for the per-wavelength flavor). When associated with a
    form that projects the luminosity along a path, the probe outputs a monochromatic luminosity
    surface density divided by the area of the unit sphere (with resulting SI units of W/m/m2/sr
    for the per-wavelength flavor). In case of parallel projection towards a distant observer, this
    quantity is equivalent to surface brightness. In other words, this probe can be used to produce
    transparent noise-free images of the dust emission. */
class SecondaryDustLuminosityProbe : public SpatialGridFormProbe
{
    /** The enumeration type indicating when probing occurs. */
    ENUM_DEF(ProbeAfter, Run, Secondary)
        ENUM_VAL(ProbeAfter, Run, "after the complete simulation run")
        ENUM_VAL(ProbeAfter, Secondary, "after each iteration over secondary emission")
    ENUM_END()

    ITEM_CONCRETE(SecondaryDustLuminosityProbe, SpatialGridFormProbe,
                  "internal spatial grid: secondary dust luminosity")
        ATTRIBUTE_TYPE_DISPLAYED_IF(SecondaryDustLuminosityProbe, "DustMix&DustEmission")

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
