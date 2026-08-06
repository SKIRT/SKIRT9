/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "CylinderSpatialGrid.hpp"
#include "FatalError.hpp"

//////////////////////////////////////////////////////////////////////

void CylinderSpatialGrid::setupSelfBefore()
{
    SpatialGrid::setupSelfBefore();
    if (_maxRadius <= _minRadius) throw FATALERROR("The outer radius must be larger than the inner radius");
    if (_maxZ <= _minZ) throw FATALERROR("The height of the cylinder must be positive");
}

//////////////////////////////////////////////////////////////////////

Box CylinderSpatialGrid::boundingBox() const
{
    return Box(-_maxRadius, -_maxRadius, _minZ, _maxRadius, _maxRadius, _maxZ);
}

//////////////////////////////////////////////////////////////////////
