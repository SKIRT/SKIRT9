/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef SPHERESPATIALGRID_HPP
#define SPHERESPATIALGRID_HPP

#include "Array.hpp"
#include "SpatialGrid.hpp"
class Mesh;

////////////////////////////////////////////////////////////////////

/** The SphereSpatialGrid class is an abstract subclass of the general SpatialGrid class, and
    represents any spatial grid defined within a spherical configuration space, centered on the
    origin of the system. */
class SphereSpatialGrid : public SpatialGrid
{
    ITEM_ABSTRACT(SphereSpatialGrid, SpatialGrid, "a spatial grid bounded by a sphere")

        PROPERTY_DOUBLE(minRadius, "the inner radius of the grid")
        ATTRIBUTE_QUANTITY(minRadius, "length")
        ATTRIBUTE_MIN_VALUE(minRadius, "[0")
        ATTRIBUTE_DEFAULT_VALUE(minRadius, "0")

        PROPERTY_DOUBLE(maxRadius, "the outer radius of the grid")
        ATTRIBUTE_QUANTITY(maxRadius, "length")
        ATTRIBUTE_MIN_VALUE(maxRadius, "]0")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function verifies the characteristics of the grid. */
    void setupSelfBefore() override;

    //======================== Other Functions =======================

public:
    /** This function returns the bounding box that encloses the grid. */
    Box boundingBox() const override;

    //======================== Helper Functions =======================

protected:
    /** This function returns \f$x^3\f$. It is intended as a convenience for subclasses. */
    static double pow3(double x) { return x * x * x; }

    /** This function returns \f$x_1^3 - x_0^3 = (x_1-x_0)(x_1^2 + x_1 x_0 + x_0^2)\f$. It is
        intended as a convenience for subclasses. */
    static double pow3(double x0, double x1) { return (x1 - x0) * (x1 * x1 + x1 * x0 + x0 * x0); }

    /** This function initializes a polar-inclination grid from the specified mesh, making sure
        that there is exactly one boundary corresponding to the equatorial plane
        (\f$\theta=\pi/2\f$ or equivalently \f$\cos\theta=0\f$). If the mesh does not include such
        as a boundary, one is automatically added. If the mesh includes multiple boundaries very
        close to the equatorial plane, a fatal error is thrown.

        The function stores the resulting boundaries and the corresponding cosine values in the
        provided arrays, resizing them appropriately, and returns the number of boundaries. */
    static int initPolarGrid(Array& thetav, Array& cosv, const Mesh* mesh);
};

////////////////////////////////////////////////////////////////////

#endif
