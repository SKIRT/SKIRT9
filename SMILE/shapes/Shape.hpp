/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef SHAPE_HPP
#define SHAPE_HPP

#include "ShapeItem.hpp"

////////////////////////////////////////////////////////////////////

/** The Shape class is the base class for all shapes that may occur in a shape hierarchy, including
    groups, decorators and actual shapes, but excluding the ShapeCanvas class (an instance of which
    sits at the root of a complete shape hierarchy). */
class Shape : public ShapeItem
{
    ITEM_ABSTRACT(Shape, ShapeItem, "a shape")
    ITEM_END()
};

////////////////////////////////////////////////////////////////////

#endif
