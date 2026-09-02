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

    // copy the configured values into our private Box
    _box = Box(_minX, _minY, _minZ, _maxX, _maxY, _maxZ);

    // verify the volume of the box
    if (_box.xwidth() <= 0) throw FATALERROR("The extent of the box should be positive in the X direction");
    if (_box.ywidth() <= 0) throw FATALERROR("The extent of the box should be positive in the Y direction");
    if (_box.zwidth() <= 0) throw FATALERROR("The extent of the box should be positive in the Z direction");

    // compute the average density
    _rho = 1. / _box.volume();
}

//////////////////////////////////////////////////////////////////////

double UniformBoxGeometry::density(Position bfr) const
{
    return _box.contains(bfr) ? _rho : 0.;
}

//////////////////////////////////////////////////////////////////////

Position UniformBoxGeometry::generatePosition() const
{
    return random()->position(_box);
}

//////////////////////////////////////////////////////////////////////

double UniformBoxGeometry::SigmaX() const
{
    if (_box.ymin() > 0. || _box.ymax() < 0. || _box.zmin() > 0. || _box.zmax() < 0.) return 0.;
    return 1. / (_box.ywidth() * _box.zwidth());
}

//////////////////////////////////////////////////////////////////////

double UniformBoxGeometry::SigmaY() const
{
    if (_box.xmin() > 0. || _box.xmax() < 0. || _box.zmin() > 0. || _box.zmax() < 0.) return 0.;
    return 1. / (_box.xwidth() * _box.zwidth());
}

//////////////////////////////////////////////////////////////////////

double UniformBoxGeometry::SigmaZ() const
{
    if (_box.xmin() > 0. || _box.xmax() < 0. || _box.ymin() > 0. || _box.ymax() < 0.) return 0.;
    return 1. / (_box.xwidth() * _box.ywidth());
}

//////////////////////////////////////////////////////////////////////
