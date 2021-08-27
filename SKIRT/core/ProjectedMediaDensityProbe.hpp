/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef PROJECTEDMEDIADENSITYPROBE_HPP
#define PROJECTEDMEDIADENSITYPROBE_HPP

#include "Probe.hpp"

////////////////////////////////////////////////////////////////////

/** ProjectedMediaDensityProbe outputs a projection of the grid-discretized media density in the
    simulation on a plane along a given line of sight. The viewing angle, field of view and pixel
    resolution are configured in the same way as for a FrameInstrument so that it is easy to obtain
    a mass density projection matching the position of an instrument.

    A separate FITS file is produced for each material type (dust, electrons, or gas), if present.
    The file for dust provides the dust surface mass density in units of dust mass per area. The
    files for electrons and gas provide the electron or hydrogen column number density in units of
    number per area. Multiple media components containing the same material type are combined,
    regardless of their ordering in the configuration. */
class ProjectedMediaDensityProbe : public Probe
{
    /** The enumeration type indicating when probing occurs. */
    ENUM_DEF(ProbeAfter, Setup, Run)
        ENUM_VAL(ProbeAfter, Setup, "after setup")
        ENUM_VAL(ProbeAfter, Run, "after the complete simulation run")
    ENUM_END()

    ITEM_CONCRETE(ProjectedMediaDensityProbe, Probe,
                  "Parallel-projection of the media density toward an arbitrary line of sight")
        ATTRIBUTE_TYPE_DISPLAYED_IF(ProjectedMediaDensityProbe, "Level2&Medium&SpatialGrid")

        PROPERTY_DOUBLE(inclination, "the inclination angle θ of the projection")
        ATTRIBUTE_QUANTITY(inclination, "posangle")
        ATTRIBUTE_MIN_VALUE(inclination, "0 deg")
        ATTRIBUTE_MAX_VALUE(inclination, "180 deg")
        ATTRIBUTE_DEFAULT_VALUE(inclination, "0 deg")
        ATTRIBUTE_DISPLAYED_IF(inclination, "Dimension2|Dimension3")

        PROPERTY_DOUBLE(azimuth, "the azimuth angle φ of the projection")
        ATTRIBUTE_QUANTITY(azimuth, "posangle")
        ATTRIBUTE_MIN_VALUE(azimuth, "-360 deg")
        ATTRIBUTE_MAX_VALUE(azimuth, "360 deg")
        ATTRIBUTE_DEFAULT_VALUE(azimuth, "0 deg")
        ATTRIBUTE_DISPLAYED_IF(azimuth, "Dimension3")

        PROPERTY_DOUBLE(roll, "the roll angle ω of the projection")
        ATTRIBUTE_QUANTITY(roll, "posangle")
        ATTRIBUTE_MIN_VALUE(roll, "-360 deg")
        ATTRIBUTE_MAX_VALUE(roll, "360 deg")
        ATTRIBUTE_DEFAULT_VALUE(roll, "0 deg")
        ATTRIBUTE_DISPLAYED_IF(roll, "Level2&(Dimension2|Dimension3)")

        PROPERTY_DOUBLE(fieldOfViewX, "the total field of view in the horizontal direction")
        ATTRIBUTE_QUANTITY(fieldOfViewX, "length")
        ATTRIBUTE_MIN_VALUE(fieldOfViewX, "]0")

        PROPERTY_INT(numPixelsX, "the number of pixels in the horizontal direction")
        ATTRIBUTE_MIN_VALUE(numPixelsX, "1")
        ATTRIBUTE_MAX_VALUE(numPixelsX, "10000")
        ATTRIBUTE_DEFAULT_VALUE(numPixelsX, "250")

        PROPERTY_DOUBLE(centerX, "the center of the frame in the horizontal direction")
        ATTRIBUTE_QUANTITY(centerX, "length")
        ATTRIBUTE_DEFAULT_VALUE(centerX, "0")
        ATTRIBUTE_DISPLAYED_IF(centerX, "Level2")

        PROPERTY_DOUBLE(fieldOfViewY, "the total field of view in the vertical direction")
        ATTRIBUTE_QUANTITY(fieldOfViewY, "length")
        ATTRIBUTE_MIN_VALUE(fieldOfViewY, "]0")

        PROPERTY_INT(numPixelsY, "the number of pixels in the vertical direction")
        ATTRIBUTE_MIN_VALUE(numPixelsY, "1")
        ATTRIBUTE_MAX_VALUE(numPixelsY, "10000")
        ATTRIBUTE_DEFAULT_VALUE(numPixelsY, "250")

        PROPERTY_DOUBLE(centerY, "the center of the frame in the vertical direction")
        ATTRIBUTE_QUANTITY(centerY, "length")
        ATTRIBUTE_DEFAULT_VALUE(centerY, "0")
        ATTRIBUTE_DISPLAYED_IF(centerY, "Level2")

        PROPERTY_ENUM(probeAfter, ProbeAfter, "when to probe the medium state")
        ATTRIBUTE_DEFAULT_VALUE(probeAfter, "Setup")
        ATTRIBUTE_DISPLAYED_IF(probeAfter, "HasDynamicState")

    ITEM_END()

    //======================== Other Functions =======================

public:
    /** This function performs probing after setup. It produces output only if the \em
        probeAfter property is set to Setup. */
    void probeSetup() override;

    /** This function performs probing after all photon packets have been emitted and detected. It
        produces output only if the \em probeAfter property is set to Run. */
    void probeRun() override;

private:
    /** This function performs the probing; it is called from probeSetup() or probeRun() depending
        on the value of the \em probeAfter property. */
    void probe();
};

////////////////////////////////////////////////////////////////////

#endif
