/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef ABSTRACTPLANARCUTSPROBE_HPP
#define ABSTRACTPLANARCUTSPROBE_HPP

#include "Probe.hpp"

////////////////////////////////////////////////////////////////////

/** AbstractPlanarCutsProbe is a base class for probes that need configuration options for planar
    cuts parallel to the coordinate planes. The provided options include the position of the planar
    cuts and the number of pixels in each spatial direction. */
class AbstractPlanarCutsProbe : public Probe
{
    ITEM_ABSTRACT(AbstractPlanarCutsProbe, Probe, "a probe using configurable planar cuts")

        PROPERTY_DOUBLE(positionX, "the x position of the cut parallel to the yz plane")
        ATTRIBUTE_QUANTITY(positionX, "length")
        ATTRIBUTE_DEFAULT_VALUE(positionX, "0")

        PROPERTY_DOUBLE(positionY, "the y position of the cut parallel to the xz plane")
        ATTRIBUTE_QUANTITY(positionY, "length")
        ATTRIBUTE_DEFAULT_VALUE(positionY, "0")

        PROPERTY_DOUBLE(positionZ, "the z position of the cut parallel to the xy plane")
        ATTRIBUTE_QUANTITY(positionZ, "length")
        ATTRIBUTE_DEFAULT_VALUE(positionZ, "0")

        PROPERTY_INT(numPixelsX, "the number of pixels in the x direction")
        ATTRIBUTE_MIN_VALUE(numPixelsX, "1")
        ATTRIBUTE_MAX_VALUE(numPixelsX, "10000")
        ATTRIBUTE_DEFAULT_VALUE(numPixelsX, "1024")

        PROPERTY_INT(numPixelsY, "the number of pixels in the y direction")
        ATTRIBUTE_MIN_VALUE(numPixelsY, "1")
        ATTRIBUTE_MAX_VALUE(numPixelsY, "10000")
        ATTRIBUTE_DEFAULT_VALUE(numPixelsY, "1024")

        PROPERTY_INT(numPixelsZ, "the number of pixels in the z direction")
        ATTRIBUTE_MIN_VALUE(numPixelsZ, "1")
        ATTRIBUTE_MAX_VALUE(numPixelsZ, "10000")
        ATTRIBUTE_DEFAULT_VALUE(numPixelsZ, "1024")

    ITEM_END()
};

////////////////////////////////////////////////////////////////////

#endif
