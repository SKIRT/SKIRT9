/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef PARALLELPROJECTIONFORM_HPP
#define PARALLELPROJECTIONFORM_HPP

#include "GenericForm.hpp"

//////////////////////////////////////////////////////////////////////

/** ParallelProjectionForm represents a generic probe form. Refer to the ProbeFormBridge class for
    more information about probes and forms.

    This particular form outputs a FITS file representing a projection of the quantity being probed
    on a plane along a given line of sight. The projection plane is considered to lie well outside
    of the model. The viewing angle, field of view and pixel resolution are configured in the same
    way as for a FrameInstrument so that it is easy to obtain a projection matching the position of
    an instrument.

    The FITS file contains a number of image frames corresponding to the number of components in a
    value of the quantity being probed (i.e. 1 for scalar quantities, 3 for vector quantities, and
    N for compound quantities). In case of a vector quantity, the three image frames represent
    the velocity vector components in the frame of the projection, i.e. the two components
    projected on the x and y axes of the projection plane and the component perpendicular to it,
    where positive values indicate vectors pointing away from the viewer. */
class ParallelProjectionForm : public GenericForm
{
    ITEM_CONCRETE(ParallelProjectionForm, GenericForm, "parallel projection on a distant plane")
        ATTRIBUTE_TYPE_DISPLAYED_IF(ParallelProjectionForm, "Level2")

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

    ITEM_END()

public:
    /** This function causes the form to output file(s) as described in the class header for the
        quantity being probed according to the information provided by the specified
        ProbeFormBridge instance. */
    void writeQuantity(const ProbeFormBridge* bridge) const override;
};

//////////////////////////////////////////////////////////////////////

#endif
