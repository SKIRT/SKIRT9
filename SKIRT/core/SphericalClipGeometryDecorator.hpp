/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef SPHERICALCLIPGEOMETRYDECORATOR_HPP
#define SPHERICALCLIPGEOMETRYDECORATOR_HPP

#include "ClipGeometryDecorator.hpp"

////////////////////////////////////////////////////////////////////

/** The SphericalClipGeometryDecorator class is a decorator that adjusts another geometry by
    setting the density equal to zero inside or outside a sphere with given position and radius.
    The dimension of the geometry implemented by a SphericalClipGeometryDecorator object depends on
    the symmetries of the geometry being decorated and on the position of the clipping sphere. */
class SphericalClipGeometryDecorator : public ClipGeometryDecorator
{
    ITEM_CONCRETE(SphericalClipGeometryDecorator, ClipGeometryDecorator,
                  "a decorator that clips another geometry using a sphere")

        PROPERTY_DOUBLE(clipRadius, "the radius of the clipping sphere")
        ATTRIBUTE_QUANTITY(clipRadius, "length")
        ATTRIBUTE_MIN_VALUE(clipRadius, "[0")

        PROPERTY_DOUBLE(centerX, "the x coordinate of the sphere's center")
        ATTRIBUTE_QUANTITY(centerX, "length")
        ATTRIBUTE_DEFAULT_VALUE(centerX, "0")
        ATTRIBUTE_DISPLAYED_IF(centerX, "Level2")
        ATTRIBUTE_INSERT(centerX, "centerX:Dimension3")

        PROPERTY_DOUBLE(centerY, "the y coordinate of the sphere's center")
        ATTRIBUTE_QUANTITY(centerY, "length")
        ATTRIBUTE_DEFAULT_VALUE(centerY, "0")
        ATTRIBUTE_DISPLAYED_IF(centerY, "Level2")
        ATTRIBUTE_INSERT(centerY, "centerY:Dimension3")

        PROPERTY_DOUBLE(centerZ, "the z coordinate of the sphere's center")
        ATTRIBUTE_QUANTITY(centerZ, "length")
        ATTRIBUTE_DEFAULT_VALUE(centerZ, "0")
        ATTRIBUTE_DISPLAYED_IF(centerZ, "Level2")
        ATTRIBUTE_INSERT(centerZ, "centerZ:Dimension2")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function calculates some frequently used values. */
    void setupSelfBefore() override;

    //======================== Other Functions =======================

public:
    /** This function returns the dimension of the geometry, which is the larger of two dimensions:
        the dimension of the geometry being decorated and the dimension of the clipping sphere. The
        dimension of the clipping sphere is 1 if its center is at the origin, 2 if the center is on
        the Z-axis, and 3 if the center is elsewhere. */
    int dimension() const override;

    /** This function returns true if the specified position is inside the sphere defined by the
        properties of this class. */
    bool inside(Position bfr) const override;

    //======================== Data Members ========================

private:
    // values initialized during setup
    Position _center;
    double _clipRadiusSquared{0.};
};

////////////////////////////////////////////////////////////////////

#endif
