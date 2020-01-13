/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "MeshSource.hpp"
#include "FatalError.hpp"

////////////////////////////////////////////////////////////////////

void MeshSource::setupSelfBefore()
{
    ImportedSource::setupSelfBefore();

    // copy the configured values into our private Box
    _domain = Box(_minX, _minY, _minZ, _maxX, _maxY, _maxZ);

    // verify the volume of the domain
    if (_domain.xwidth() <= 0) throw FATALERROR("The extent of the domain should be positive in the X direction");
    if (_domain.ywidth() <= 0) throw FATALERROR("The extent of the domain should be positive in the Y direction");
    if (_domain.zwidth() <= 0) throw FATALERROR("The extent of the domain should be positive in the Z direction");
}

////////////////////////////////////////////////////////////////////
