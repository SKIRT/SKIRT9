/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef MERIODIONALDUSTTEMPERATURECUTPROBE_HPP
#define MERIODIONALDUSTTEMPERATURECUTPROBE_HPP

#include "Probe.hpp"

////////////////////////////////////////////////////////////////////

/** MeridionalDustTemperatureCutProbe outputs a column text file (named
    <tt>prefix_probe_T.dat</tt>) listing an indicative dust temperature along a meridian line
    (i.e., one half of a great circle) through the medium. The half-circle is centered on the
    origin of the model coordinate system, lies in a meridional plane, and spans all inclinations
    (\f$0\leq\theta\leq\pi\f$). The radius \f$r\f$ and the azimuth \f$\varphi\f$ of the meridian
    can be configured, and the number of equidistant samples to be taken along the meridian can be
    specified as well. The resolution is limited only by the resolution of the simulation's spatial
    grid.

    The output file contains a line for each sample along the meridian. The first column specifies
    the inclination \f$\theta\f$, and the second column lists the corresponding indicative dust
    temperature.

    The indicative dust temperature for a particular sample position is obtained as follows. First
    the probe determines the cell in the simulation's spatial grid containing the position. For
    each material mix of type dust present in the cell, or if applicable, for each dust population
    in these mixes, the probe calculates the equilibrium temperature that would be reached when the
    dust is embedded in the radiation field tracked by the simulation for the cell. This is
    achieved by solving the energy balance equation under LTE (local thermal equilibrium)
    assumptions. The resulting temperatures are finally averaged over the dust populations in each
    mix (weighed by the relative mass in the mix) and over all dust components present in the
    spatial cell (weighed by relative mass in the cell).

    Note that the indicative dust temperature does not correspond to a physical temperature. The
    LTE assumption is almost certainly unjustified for a relevant portion of the dust grains
    (depending on the embedding radiation field), and even when ignoring this problem, averaging
    temperatures over dust populations and dust mixes has no clear-cut physical interpretation. */
class MeridionalDustTemperatureCutProbe : public Probe
{
    ITEM_CONCRETE(MeridionalDustTemperatureCutProbe, Probe, "the indicative dust temperature along a meridian line")
        ATTRIBUTE_TYPE_DISPLAYED_IF(MeridionalDustTemperatureCutProbe, "Level2&Dust&SpatialGrid&RadiationField")

    PROPERTY_INT(numSamples, "the number of samples along the meridian")
        ATTRIBUTE_MIN_VALUE(numSamples, "3")
        ATTRIBUTE_MAX_VALUE(numSamples, "100000")
        ATTRIBUTE_DEFAULT_VALUE(numSamples, "250")

    PROPERTY_DOUBLE(radius, "the radius of the circle containing the meridian")
        ATTRIBUTE_QUANTITY(radius, "length")
        ATTRIBUTE_MIN_VALUE(radius, "]0")

    PROPERTY_DOUBLE(azimuth, "the azimuth angle φ of the meridian")
        ATTRIBUTE_QUANTITY(azimuth, "posangle")
        ATTRIBUTE_MIN_VALUE(azimuth, "-360 deg")
        ATTRIBUTE_MAX_VALUE(azimuth, "360 deg")
        ATTRIBUTE_DEFAULT_VALUE(azimuth, "0 deg")
        ATTRIBUTE_DISPLAYED_IF(azimuth, "Dimension3")

    ITEM_END()

    //======================== Other Functions =======================

public:
    /** This function performs probing after all photon packets have been emitted and detected. */
    void probeRun() override;
};

////////////////////////////////////////////////////////////////////

#endif
