/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef LINEARDUSTTEMPERATURECUTPROBE_HPP
#define LINEARDUSTTEMPERATURECUTPROBE_HPP

#include "Probe.hpp"

////////////////////////////////////////////////////////////////////

/** LinearDustTemperatureCutProbe outputs a column text file (named <tt>prefix_probe_T.dat</tt>)
    listing an indicative dust temperature along an arbitrary line segment throught the medium. The
    starting and ending point of the line segment can be configured using regular model
    coordinates, and the number of equidistant samples to be taken along the line can be specified
    as well. The resolution is limited only by the resolution of the simulation's spatial grid.

    The output file contains a line for each sample along the line segment. The first column
    specifies the distance along the line from the starting point, and the second column lists the
    indicative dust temperature.

    The indicative dust temperature for a particular sample position is obtained as follows. First
    the probe determines the cell in the simulation's spatial grid containing the position. The
    probe then calculates the indicative dust temperature for that cell by averaging the LTE
    equilibrium temperatures for the various dust mixes present in the cell. Note that the
    indicative dust temperature does not really correspond to a physical temperature. For more
    information about the indicative dust temperature, refer to the
    MediumSystem::indicativeDustTemperature() function. */
class LinearDustTemperatureCutProbe : public Probe
{
    ITEM_CONCRETE(LinearDustTemperatureCutProbe, Probe, "the indicative dust temperature along a given line segment")
        ATTRIBUTE_TYPE_DISPLAYED_IF(LinearDustTemperatureCutProbe,
                                    "Level2&Dust&SpatialGrid&RadiationField&Panchromatic")

    PROPERTY_INT(numSamples, "the number of samples along the line segment")
        ATTRIBUTE_MIN_VALUE(numSamples, "3")
        ATTRIBUTE_MAX_VALUE(numSamples, "100000")
        ATTRIBUTE_DEFAULT_VALUE(numSamples, "250")

    PROPERTY_DOUBLE(startX, "the position of the starting point, x component")
        ATTRIBUTE_QUANTITY(startX, "length")

    PROPERTY_DOUBLE(startY, "the position of the starting point, y component")
        ATTRIBUTE_QUANTITY(startY, "length")

    PROPERTY_DOUBLE(startZ, "the position of the starting point, z component")
        ATTRIBUTE_QUANTITY(startZ, "length")

    PROPERTY_DOUBLE(endX, "the position of the ending point, x component")
        ATTRIBUTE_QUANTITY(endX, "length")

    PROPERTY_DOUBLE(endY, "the position of the ending point, y component")
        ATTRIBUTE_QUANTITY(endY, "length")

    PROPERTY_DOUBLE(endZ, "the position of the ending point, z component")
        ATTRIBUTE_QUANTITY(endZ, "length")

    ITEM_END()

    //======================== Other Functions =======================

public:
    /** This function performs probing after all photon packets have been emitted and detected. */
    void probeRun() override;
};

////////////////////////////////////////////////////////////////////

#endif
