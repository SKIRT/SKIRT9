/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef BOXGEOMETRY_HPP
#define BOXGEOMETRY_HPP

#include "GenGeometry.hpp"
#include "Box.hpp"

////////////////////////////////////////////////////////////////////

/** The BoxGeometry class is an abstract subclass of the Geometry class, and represents any
    geometry defined within a cuboidal volume with faces that are aligned with the planes of the
    coordinate system (a box). The class also inherits from the Box class. */
class BoxGeometry : public GenGeometry, public Box
{
    ITEM_ABSTRACT(BoxGeometry, GenGeometry, "a geometry bounded by a box")

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

protected:
    /** This function verifies that the box has a positive volume. */
    void setupSelfBefore() override;
};

////////////////////////////////////////////////////////////////////

#endif
