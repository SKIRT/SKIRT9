/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "BoxGeometry.hpp"
#include "FatalError.hpp"

//////////////////////////////////////////////////////////////////////

void BoxGeometry::setupSelfBefore()
{
    GenGeometry::setupSelfBefore();

    // verify that the box is not empty
    if (_maxX <= _minX) throw FATALERROR("The extent of the box should be positive in the X direction");
    if (_maxY <= _minY) throw FATALERROR("The extent of the box should be positive in the Y direction");
    if (_maxZ <= _minZ) throw FATALERROR("The extent of the box should be positive in the Z direction");

    // copy the configured values into our inherited Box
    setExtent(_minX, _minY, _minZ, _maxX, _maxY, _maxZ);
}

//////////////////////////////////////////////////////////////////////
