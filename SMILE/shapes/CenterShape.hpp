/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef CENTERSHAPE_HPP
#define CENTERSHAPE_HPP

#include "Shape.hpp"

////////////////////////////////////////////////////////////////////

/** The CenterShape class is the base class for all low-level shapes that actually draw marks
    surrounding a central position. The CenterShape class inherits from the Shape class, just like
    all shapes should, and offers properties that specify the central position of the shape. */
class CenterShape : public Shape
{
    ITEM_ABSTRACT(CenterShape, Shape, "a shape with a central position")
        PROPERTY_DOUBLE(x, "the horizontal coordinate of the shape center")
        ATTRIBUTE_QUANTITY(x, "length")
        ATTRIBUTE_MIN_VALUE(x, "0 m")
        ATTRIBUTE_MAX_VALUE(x, "1 m")
        PROPERTY_DOUBLE(y, "the vertical coordinate of the shape center")
        ATTRIBUTE_QUANTITY(y, "length")
        ATTRIBUTE_MIN_VALUE(y, "0 m")
        ATTRIBUTE_MAX_VALUE(y, "1 m")
    ITEM_END()
};

////////////////////////////////////////////////////////////////////

#endif
