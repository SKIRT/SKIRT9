/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "BoxClipGeometryDecorator.hpp"
#include "FatalError.hpp"

////////////////////////////////////////////////////////////////////

void BoxClipGeometryDecorator::setupSelfBefore()
{
    ClipGeometryDecorator::setupSelfBefore();

    // copy the configured values into our private Box
    _box = Box(_minX, _minY, _minZ, _maxX, _maxY, _maxZ);

    // verify the volume of the box
    if (_box.xwidth() <= 0) throw FATALERROR("The extent of the box should be positive in the X direction");
    if (_box.ywidth() <= 0) throw FATALERROR("The extent of the box should be positive in the Y direction");
    if (_box.zwidth() <= 0) throw FATALERROR("The extent of the box should be positive in the Z direction");
}

////////////////////////////////////////////////////////////////////

int BoxClipGeometryDecorator::dimension() const
{
    return 3;
}

////////////////////////////////////////////////////////////////////

bool BoxClipGeometryDecorator::inside(Position bfr) const
{
    return _box.contains(bfr);
}

////////////////////////////////////////////////////////////////////
