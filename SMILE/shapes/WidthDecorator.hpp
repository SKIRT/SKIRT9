/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef WIDTHDECORATOR_HPP
#define WIDTHDECORATOR_HPP

#include "ShapeDecorator.hpp"

////////////////////////////////////////////////////////////////////

/** The WidthDecorator class represents a decorator that sets the line width for the shape or group
    of shapes being decorated. */
class WidthDecorator : public ShapeDecorator
{
    ITEM_CONCRETE(WidthDecorator, ShapeDecorator, "a decorator that sets the line width of a shape or group of shapes")
        ATTRIBUTE_TYPE_ALLOWED_IF(WidthDecorator, "!PolygonShape")
        PROPERTY_DOUBLE(width, "the line width of the decorated shape or group of shapes")
        ATTRIBUTE_MIN_VALUE(width, "0")
        ATTRIBUTE_MAX_VALUE(width, "1")
        ATTRIBUTE_DEFAULT_VALUE(width, "0.01")
    ITEM_END()

    // ================== Drawing ==================

protected:
    /** This function sets the line width in the current graphics state by calling the setWidth()
        function. */
    void paintSelf() override;
};

////////////////////////////////////////////////////////////////////

#endif
