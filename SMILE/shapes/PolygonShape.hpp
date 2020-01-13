/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef POLYGONSHAPE_HPP
#define POLYGONSHAPE_HPP

#include "CenterShape.hpp"

////////////////////////////////////////////////////////////////////

/** The PolygonShape class represents a regular polygon that draws its \f$N \leq 3\f$ sides. The
    class inherits from the CenterShape class. All corner points lie at the same radius \f$r>0\f$
    from the center of the shape. One of the corner points lies on the horizontal or vertical line
    through the center depending on the value of the \em align property. */
class PolygonShape : public CenterShape
{
    ITEM_CONCRETE(PolygonShape, CenterShape, "a regular polygon")
        PROPERTY_DOUBLE(radius, "the radius, i.e. the distance from each corner to the center")
        ATTRIBUTE_QUANTITY(radius, "length")
        ATTRIBUTE_MIN_VALUE(radius, "0 m")
        ATTRIBUTE_MAX_VALUE(radius, "1 m")
        PROPERTY_INT(numSides, "the number of sides (or corners)")
        ATTRIBUTE_MIN_VALUE(numSides, "3")
        ATTRIBUTE_MAX_VALUE(numSides, "99")
        PROPERTY_BOOL(align, "align the first corner with the vertical line through the center")
        ATTRIBUTE_DEFAULT_VALUE(align, "false")
    ITEM_END()

    // ================== Drawing ==================

protected:
    /** This function paints the appropriate marks for the regular polygon shape by calling the
        drawLine() function for each side. */
    void paintSelf() override;
};

////////////////////////////////////////////////////////////////////

#endif
