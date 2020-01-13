/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "SphericalClipGeometryDecorator.hpp"

////////////////////////////////////////////////////////////////////

void SphericalClipGeometryDecorator::setupSelfBefore()
{
    ClipGeometryDecorator::setupSelfBefore();

    // calculate some frequently used values
    _center = Position(_centerX, _centerY, _centerZ);
    _clipRadiusSquared = _clipRadius * _clipRadius;
}

////////////////////////////////////////////////////////////////////

int SphericalClipGeometryDecorator::dimension() const
{
    int clipDimension = 1;
    if (_centerZ) clipDimension = 2;
    if (_centerX || _centerY) clipDimension = 3;
    return std::max(geometry()->dimension(), clipDimension);
}

////////////////////////////////////////////////////////////////////

bool SphericalClipGeometryDecorator::inside(Position bfr) const
{
    return (bfr - _center).norm2() <= _clipRadiusSquared;
}

////////////////////////////////////////////////////////////////////
