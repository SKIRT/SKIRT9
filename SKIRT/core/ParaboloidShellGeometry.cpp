/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "ParaboloidShellGeometry.hpp"
#include "FatalError.hpp"
#include "Random.hpp"

//////////////////////////////////////////////////////////////////////

void ParaboloidShellGeometry::setupSelfBefore()
{
    AxGeometry::setupSelfBefore();

    // verify property values
    if (_DeltaOut <= _DeltaIn) throw FATALERROR("the outer opening angle should be larger than the inner");
    if (_zOut >= _zIn) throw FATALERROR("the outer paraboloid vertical offset should be smaller than the inner");

    // determine the inner paraboloid height
    _zmaxIn = _Din * cos(_DeltaIn);

    // determine the inner paraboloid radius
    _Rin = _Din * sin(_DeltaIn);

    // determine the curvature level of the inner paraboloid wall
    _aIn = sqrt(tan(_DeltaIn) * _Rin);

    // determine the outer paraboloid height
    _zmaxOut = _zmaxIn + _zIn - _zOut;

    // determine the outer paraboloid radius
    _Rout = _zmaxOut * tan(_DeltaOut);

    // determine the curvature level of the outer paraboloid wall
    _aOut = sqrt(tan(_DeltaOut) * _Rout);

    // determine the normalization factor
    _A = 1. / M_PI / (_Rout * _Rout * _zmaxOut - _Rin * _Rin * _zmaxIn);
}

////////////////////////////////////////////////////////////////////

double ParaboloidShellGeometry::density(double R, double z) const
{
    if (z > (_zmaxIn + _zIn)) return 0.;
    if (z < -(_zmaxIn + _zIn)) return 0.;
    if ((z >= 0.) && (z < ((R * R) / (_aOut * _aOut) + _zOut))) return 0.;
    if ((z >= 0.) && (z > ((R * R) / (_aIn * _aIn) + _zIn))) return 0.;
    if ((z < 0.) && (z > -((R * R) / (_aOut * _aOut) + _zOut))) return 0.;
    if ((z < 0.) && (z < -((R * R) / (_aIn * _aIn) + _zIn))) return 0.;
    return _A;
}

//////////////////////////////////////////////////////////////////////

Position ParaboloidShellGeometry::generatePosition() const
{
    while (true)
    {
        double R = _Rout * sqrt(random()->uniform());
        double z = (_zmaxIn + _zIn) * (2. * random()->uniform() - 1.);
        if (density(R, z))
        {
            double phi = 2. * M_PI * random()->uniform();
            return Position(R, phi, z, Position::CoordinateSystem::CYLINDRICAL);
        }
    }
}

//////////////////////////////////////////////////////////////////////

double ParaboloidShellGeometry::SigmaR() const
{
    return 0.;
}

//////////////////////////////////////////////////////////////////////

double ParaboloidShellGeometry::SigmaZ() const
{
    return 2. * _A * (_zIn - _zOut);
}

//////////////////////////////////////////////////////////////////////
