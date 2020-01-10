/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef RECTANGLESHAPE_HPP
#define RECTANGLESHAPE_HPP

#include "CenterShape.hpp"

////////////////////////////////////////////////////////////////////

/** The RectangleShape class represents a rectangular shape (lined up with the coordinate axes)
    that draws its four sides. The RectangleShape class inherits from the CenterShape class, and
    offers additional properties that specify the size of the rectangle. */
class RectangleShape : public CenterShape
{
    ITEM_CONCRETE(RectangleShape, CenterShape, "a rectangular shape")
        PROPERTY_DOUBLE(width, "the width of the rectangle")
        ATTRIBUTE_QUANTITY(width, "length")
        ATTRIBUTE_MIN_VALUE(width, "0 m")
        ATTRIBUTE_MAX_VALUE(width, "1 m")
        PROPERTY_DOUBLE(height, "the height of the rectangle")
        ATTRIBUTE_QUANTITY(height, "length")
        ATTRIBUTE_MIN_VALUE(height, "0 m")
        ATTRIBUTE_MAX_VALUE(height, "1 m")
    ITEM_END()

    // ================== Drawing ==================

protected:
    /** This function paints the appropriate marks for the rectangular shape by calling the
        drawLine() function for each side. */
    void paintSelf() override;
};

////////////////////////////////////////////////////////////////////

#endif
