/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef SHAPEDECORATOR_HPP
#define SHAPEDECORATOR_HPP

#include "Shape.hpp"

////////////////////////////////////////////////////////////////////

/** ShapeDecorator is the base class for all shape decorators, i.e. shapes that adjust another
    shape (or group of shapes). The specific form of adjustment is defined by each ShapeDecorator
    subclass. */
class ShapeDecorator : public Shape
{
    ITEM_ABSTRACT(ShapeDecorator, Shape, "a shape decorator adjusting a shape or group of shapes")
        PROPERTY_ITEM(shape, Shape, "the shape or group of shapes being decorated")
    ITEM_END()
};

////////////////////////////////////////////////////////////////////

#endif
