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
    the probe determines the cell in the simulation's spatial grid containing the position. The
    probe then calculates the indicative dust temperature for that cell by averaging the LTE
    equilibrium temperatures for the various dust mixes present in the cell. Note that the
    indicative dust temperature does not really correspond to a physical temperature. For more
    information about the indicative dust temperature, refer to the
    MediumSystem::indicativeDustTemperature() function. */
class MeridionalDustTemperatureCutProbe : public Probe
{
    ITEM_CONCRETE(MeridionalDustTemperatureCutProbe, Probe, "the indicative dust temperature along a meridian line")
        ATTRIBUTE_TYPE_DISPLAYED_IF(MeridionalDustTemperatureCutProbe,
                                    "Level2&Dust&SpatialGrid&RadiationField&Panchromatic")

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
