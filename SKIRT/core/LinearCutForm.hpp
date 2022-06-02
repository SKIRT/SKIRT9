/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef LINEARCUTFORM_HPP
#define LINEARCUTFORM_HPP

#include "GenericForm.hpp"

//////////////////////////////////////////////////////////////////////

/** LinearCutForm represents a generic probe form. Refer to the ProbeFormBridge class for more
    information about probes and forms.

    This particular form outputs a text column file listing the quantity being probed at a number
    of positions along an arbitrary line segment through the spatial domain of the simulation. The
    starting and ending point of the line segment can be configured using regular model
    coordinates, and the number of equidistant samples to be taken along the line can be specified
    as well.

    The output file contains a line for each sample along the line segment. The first column
    specifies the distance along the line from the starting point, and subsequent column(s) list
    the quantity being probed. */
class LinearCutForm : public GenericForm
{
    ITEM_CONCRETE(LinearCutForm, GenericForm, "a text column file with values along a given line segment")
        ATTRIBUTE_TYPE_DISPLAYED_IF(LinearCutForm, "Level2")

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

public:
    /** This function causes the form to output file(s) as described in the class header for the
        quantity being probed according to the information provided by the specified
        ProbeFormBridge instance. */
    void writeQuantity(const ProbeFormBridge* bridge) const override;
};

//////////////////////////////////////////////////////////////////////

#endif
