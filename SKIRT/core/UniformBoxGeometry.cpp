/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "UniformBoxGeometry.hpp"
#include "FatalError.hpp"
#include "Random.hpp"

////////////////////////////////////////////////////////////////////

void UniformBoxGeometry::setupSelfBefore()
{
    GenGeometry::setupSelfBefore();

    // verify that the box is not empty
    if (_maxX <= _minX) throw FATALERROR("The extent of the box should be positive in the X direction");
    if (_maxY <= _minY) throw FATALERROR("The extent of the box should be positive in the Y direction");
    if (_maxZ <= _minZ) throw FATALERROR("The extent of the box should be positive in the Z direction");

    // copy the configured values into our inherited Box
    setExtent(_minX, _minY, _minZ, _maxX, _maxY, _maxZ);

    // compute the average density
    _rho = 1. / volume();
}

//////////////////////////////////////////////////////////////////////

double UniformBoxGeometry::density(Position bfr) const
{
    return contains(bfr) ? _rho : 0.;
}

//////////////////////////////////////////////////////////////////////

Position UniformBoxGeometry::generatePosition() const
{
    return random()->position(extent());
}

//////////////////////////////////////////////////////////////////////

double UniformBoxGeometry::SigmaX() const
{
    return 1. / (ywidth() * zwidth());
}

//////////////////////////////////////////////////////////////////////

double UniformBoxGeometry::SigmaY() const
{
    return 1. / (xwidth() * zwidth());
}

//////////////////////////////////////////////////////////////////////

double UniformBoxGeometry::SigmaZ() const
{
    return 1. / (xwidth() * ywidth());
}

//////////////////////////////////////////////////////////////////////
