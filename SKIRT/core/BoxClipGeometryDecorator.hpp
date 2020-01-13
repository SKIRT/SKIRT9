/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef BOXCLIPGEOMETRYDECORATOR_HPP
#define BOXCLIPGEOMETRYDECORATOR_HPP

#include "Box.hpp"
#include "ClipGeometryDecorator.hpp"

////////////////////////////////////////////////////////////////////

/** The BoxClipGeometryDecorator class is a decorator that adjusts another geometry by setting the
    density equal to zero inside or outside a given cuboidal bounding box. */
class BoxClipGeometryDecorator : public ClipGeometryDecorator
{
    ITEM_CONCRETE(BoxClipGeometryDecorator, ClipGeometryDecorator,
                  "a decorator that clips another geometry using a cubodial box")
        ATTRIBUTE_TYPE_INSERT(BoxClipGeometryDecorator, "Dimension3")

        PROPERTY_DOUBLE(minX, "the start point of the box in the X direction")
        ATTRIBUTE_QUANTITY(minX, "length")

        PROPERTY_DOUBLE(maxX, "the end point of the box in the X direction")
        ATTRIBUTE_QUANTITY(maxX, "length")

        PROPERTY_DOUBLE(minY, "the start point of the box in the Y direction")
        ATTRIBUTE_QUANTITY(minY, "length")

        PROPERTY_DOUBLE(maxY, "the end point of the box in the Y direction")
        ATTRIBUTE_QUANTITY(maxY, "length")

        PROPERTY_DOUBLE(minZ, "the start point of the box in the Z direction")
        ATTRIBUTE_QUANTITY(minZ, "length")

        PROPERTY_DOUBLE(maxZ, "the end point of the box in the Z direction")
        ATTRIBUTE_QUANTITY(maxZ, "length")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

public:
    /** This function verifies that the box has a positive volume. */
    void setupSelfBefore() override;

    //======================== Other Functions =======================

public:
    /** This function returns the dimension of the geometry, which is 3 for this
        class since a box has no specific symmetry. */
    int dimension() const override;

    /** This function returns true if the specified position is inside the box defined by the
        properties of this class. */
    bool inside(Position bfr) const override;

    //======================== Data Members ========================

private:
    Box _box;
};

////////////////////////////////////////////////////////////////////

#endif
