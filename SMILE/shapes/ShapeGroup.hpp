/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef SHAPEGROUP_HPP
#define SHAPEGROUP_HPP

#include "Shape.hpp"

////////////////////////////////////////////////////////////////////

/** The ShapeGroup class is a Shape subclass that groups an arbitrary number of shapes into a
    single compound shape. The class performs pure aggregation, i.e. it does not alter or augment
    the grouped shapes in any way. */
class ShapeGroup : public Shape
{
    ITEM_CONCRETE(ShapeGroup, Shape, "a compound shape that groups other shapes")
        PROPERTY_ITEM_LIST(shapes, Shape, "the shapes in this group")
    ITEM_END()
};

////////////////////////////////////////////////////////////////////

#endif
