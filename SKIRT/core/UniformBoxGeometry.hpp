/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef UNIFORMBOXGEOMETRY_HPP
#define UNIFORMBOXGEOMETRY_HPP

#include "Box.hpp"
#include "GenGeometry.hpp"

////////////////////////////////////////////////////////////////////

/** The UniformBoxGeometry class describes a 3D geometry consisting of a uniform cuboid aligned
    with the coordinate system. */
class UniformBoxGeometry : public GenGeometry
{
    ITEM_CONCRETE(UniformBoxGeometry, GenGeometry, "a geometry with uniform density inside a box")

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
    /** This function verifies that the box has a positive volume, and calculates the (constant)
        density at each point, assuming a total mass of 1. */
    void setupSelfBefore() override;

    //======================== Other Functions =======================

public:
    /** This function returns the density \f$\rho({\bf{r}})\f$ for this geometry at the position
        \f${\bf{r}}\f$. It returns the constant value
        \f$\rho=1/(x_\text{width}\,y_\text{width}\,z_\text{width})\f$ for all points inside the
        box, and zero outside the box. */
    double density(Position bfr) const override;

    /** This function generates a random position from the geometry, by drawing a uniform random
        point inside the box. */
    Position generatePosition() const override;

    /** This function returns the X-axis surface density of the geometry, defined as the
        integration of the density along the entire X-axis. */
    double SigmaX() const override;

    /** This function returns the Y-axis surface density of the geometry, defined as the
        integration of the density along the entire Y-axis. */
    double SigmaY() const override;

    /** This function returns the Z-axis surface density of the geometry, defined as the
        integration of the density along the entire Z-axis. */
    double SigmaZ() const override;

    //======================== Data Members ========================

private:
    // data members initialized during setup
    double _rho{0.};
    Box _box;
};

////////////////////////////////////////////////////////////////////

#endif
