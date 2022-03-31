/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef ALLSKYPROJECTIONFORM_HPP
#define ALLSKYPROJECTIONFORM_HPP

#include "AllSkyProjection.hpp"
#include "GenericForm.hpp"

//////////////////////////////////////////////////////////////////////

/** AllSkyProjectionForm represents a generic probe form. Refer to the ProbeFormBridge class for
    more information about probes and forms.

    This particular form outputs a FITS file representing an all-sky projection of the quantity
    being probed as seen from the given position and centered towards the given direction. The
    position can be inside the model. The map has a 2:1 aspect ratio (with square pixels) and uses
    the user-selected projection to project the complete sky on the rectangular image. Pixels
    outside of the mapped region are set to zero.

    The FITS file contains a number of image frames corresponding to the number of components in a
    value of the quantity being probed (i.e. 1 for scalar quantities, 3 for vector quantities, and
    N for compound quantities). In case of a vector quantity, the three image frames represent the
    velocity vector components in the frame of the projection, i.e. the two components projected on
    the x and y axes of the projection plane and the component perpendicular to it, where positive
    values indicate vectors pointing away from the viewer. */
class AllSkyProjectionForm : public GenericForm
{
    ITEM_CONCRETE(AllSkyProjectionForm, GenericForm, "all-sky projection at a position inside the model")
        ATTRIBUTE_TYPE_DISPLAYED_IF(AllSkyProjectionForm, "Level2")

        PROPERTY_ITEM(projection, AllSkyProjection, "the projection used for mapping the sky to a rectangle")
        ATTRIBUTE_DEFAULT_VALUE(projection, "HammerAitoffProjection")

        PROPERTY_INT(numPixelsY, "the number of image pixels in the vertical (shortest) direction")
        ATTRIBUTE_MIN_VALUE(numPixelsY, "25")
        ATTRIBUTE_MAX_VALUE(numPixelsY, "10000")
        ATTRIBUTE_DEFAULT_VALUE(numPixelsY, "250")

        PROPERTY_DOUBLE(observerX, "the position of the observer, x component")
        ATTRIBUTE_QUANTITY(observerX, "length")

        PROPERTY_DOUBLE(observerY, "the position of the observer, y component")
        ATTRIBUTE_QUANTITY(observerY, "length")

        PROPERTY_DOUBLE(observerZ, "the position of the observer, z component")
        ATTRIBUTE_QUANTITY(observerZ, "length")

        PROPERTY_DOUBLE(crossX, "the position of the crosshair, x component")
        ATTRIBUTE_QUANTITY(crossX, "length")
        ATTRIBUTE_DEFAULT_VALUE(crossX, "0")

        PROPERTY_DOUBLE(crossY, "the position of the crosshair, y component")
        ATTRIBUTE_QUANTITY(crossY, "length")
        ATTRIBUTE_DEFAULT_VALUE(crossY, "0")

        PROPERTY_DOUBLE(crossZ, "the position of the crosshair, z component")
        ATTRIBUTE_QUANTITY(crossZ, "length")
        ATTRIBUTE_DEFAULT_VALUE(crossZ, "0")

        PROPERTY_DOUBLE(upX, "the upwards direction, x component")
        ATTRIBUTE_QUANTITY(upX, "length")
        ATTRIBUTE_DEFAULT_VALUE(upX, "0")

        PROPERTY_DOUBLE(upY, "the upwards direction, y component")
        ATTRIBUTE_QUANTITY(upY, "length")
        ATTRIBUTE_DEFAULT_VALUE(upY, "0")

        PROPERTY_DOUBLE(upZ, "the upwards direction, z component")
        ATTRIBUTE_QUANTITY(upZ, "length")
        ATTRIBUTE_DEFAULT_VALUE(upZ, "1")

    ITEM_END()

protected:
    /** This function verifies that all attribute values have been appropriately set. */
    void setupSelfBefore() override;

public:
    /** This function causes the form to output file(s) as described in the class header for the
        quantity being probed according to the information provided by the specified
        ProbeFormBridge instance. */
    void writeQuantity(const ProbeFormBridge* bridge) const override;

private:
    // aliases to discoverable data members for ease of notation
    const double& _Ox{_observerX};
    const double& _Oy{_observerY};
    const double& _Oz{_observerZ};
    const double& _Cx{_crossX};
    const double& _Cy{_crossY};
    const double& _Cz{_crossZ};
    const double& _Ux{_upX};
    const double& _Uy{_upY};
    const double& _Uz{_upZ};
};

//////////////////////////////////////////////////////////////////////

#endif
