/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "ParaboloidGeometry.hpp"
#include "FatalError.hpp"
#include "Random.hpp"

//////////////////////////////////////////////////////////////////////

void ParaboloidGeometry::setupSelfBefore()
{
    AxGeometry::setupSelfBefore();

    // determine the paraboloid radius
    _Rp = _D * sin(_Delta);

    // determine the paraboloid height
    _zmax = _D * cos(_Delta);

    // determine the curvature level of the paraboloid
    _a = sqrt(tan(_Delta) * _Rp);

    // determine the normalization factor
    _A = 1. / M_PI / (_Rp * _Rp) / _zmax;
}

////////////////////////////////////////////////////////////////////

double ParaboloidGeometry::density(double R, double z) const
{
    if (z > (_zmax + _z0)) return 0.;
    if (z < -(_zmax + _z0)) return 0.;
    if ((z >= 0.) && (z < ((R * R) / (_a * _a) + _z0))) return 0.;
    if ((z < 0.) && (z > -((R * R) / (_a * _a) + _z0))) return 0.;
    return _A;
}

//////////////////////////////////////////////////////////////////////

Position ParaboloidGeometry::generatePosition() const
{
    while (true)
    {
        double R = _Rp * sqrt(random()->uniform());
        double z = (_zmax + _z0) * (2. * random()->uniform() - 1.);
        if (density(R, z))
        {
            double phi = 2. * M_PI * random()->uniform();
            return Position(R, phi, z, Position::CoordinateSystem::CYLINDRICAL);
        }
    }
}

//////////////////////////////////////////////////////////////////////

double ParaboloidGeometry::SigmaR() const
{
    return 0.;
}

//////////////////////////////////////////////////////////////////////

double ParaboloidGeometry::SigmaZ() const
{
    return 2. * _A * _zmax;
}

//////////////////////////////////////////////////////////////////////
