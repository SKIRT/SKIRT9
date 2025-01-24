/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef SPHERESPATIALGRID_HPP
#define SPHERESPATIALGRID_HPP

#include "SpatialGrid.hpp"

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
};

////////////////////////////////////////////////////////////////////

#endif
