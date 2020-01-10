/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef CYLINDERSPATIALGRID_HPP
#define CYLINDERSPATIALGRID_HPP

#include "SpatialGrid.hpp"

////////////////////////////////////////////////////////////////////

/** The CylinderSpatialGrid class is an abstract subclass of the general SpatialGrid class, and
    represents any spatial grid defined within a cylindrical configuration space, with the symmetry
    axis the Z-axis of the system. */
class CylinderSpatialGrid : public SpatialGrid
{
    ITEM_ABSTRACT(CylinderSpatialGrid, SpatialGrid, "a spatial grid bounded by a cylinder")

        PROPERTY_DOUBLE(maxRadius, "the cylindrical radius of the grid")
        ATTRIBUTE_QUANTITY(maxRadius, "length")
        ATTRIBUTE_MIN_VALUE(maxRadius, "]0")

        PROPERTY_DOUBLE(minZ, "the start point of the cylinder in the Z direction")
        ATTRIBUTE_QUANTITY(minZ, "length")

        PROPERTY_DOUBLE(maxZ, "the end point of the cylinder in the Z direction")
        ATTRIBUTE_QUANTITY(maxZ, "length")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function verifies the characteristics of the grid. */
    void setupSelfBefore() override;

    //======================== Other Functions =======================

    /** This function returns the bounding box that encloses the grid. */
    Box boundingBox() const override;
};

////////////////////////////////////////////////////////////////////

#endif
