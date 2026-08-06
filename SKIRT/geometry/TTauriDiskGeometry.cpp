/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "TTauriDiskGeometry.hpp"
#include "FatalError.hpp"
#include "Random.hpp"
#include "SpecialFunctions.hpp"

////////////////////////////////////////////////////////////////////

void TTauriDiskGeometry::setupSelfBefore()
{
    AxGeometry::setupSelfBefore();

    // verify property values
    if (_Rout <= _Rinn) throw FATALERROR("the outer radius of the disk must be larger than the inner radius");

    // calculate cached values
    _a178 = _a - 17. / 8.;
    _glnInn = SpecialFunctions::gln(_a178, _Rinn / _Rd);
    _glnInnOut = SpecialFunctions::gln2(_a178, _Rout / _Rd, _Rinn / _Rd);
    _rho0 = 1. / (2. * pow(M_PI, 1.5) * pow(_b, -0.5) * (_Rd * _Rd * _zd) * _glnInnOut);
    _s0 = _zd / sqrt(2. * _b);
}

//////////////////////////////////////////////////////////////////////

double TTauriDiskGeometry::density(double R, double z) const
{
    if (R < _Rinn || R > _Rout) return 0.;
    double x = (z / _zd) * pow(R / _Rd, -9. / 8.);
    return _rho0 * pow(R / _Rd, -_a) * exp(-_b * x * x);
}

//////////////////////////////////////////////////////////////////////

Position TTauriDiskGeometry::generatePosition() const
{
    double R = _Rd * SpecialFunctions::gexp(_a178, _glnInn + random()->uniform() * _glnInnOut);
    double sigma = _s0 * pow(R / _Rd, 9. / 8.);
    double z = random()->gauss() * sigma;
    double phi = 2. * M_PI * random()->uniform();
    return Position(R, phi, z, Position::CoordinateSystem::CYLINDRICAL);
}

//////////////////////////////////////////////////////////////////////

double TTauriDiskGeometry::SigmaR() const
{
    return _rho0 * _Rd * SpecialFunctions::gln2(_a, _Rout / _Rd, _Rinn / _Rd);
}

//////////////////////////////////////////////////////////////////////

double TTauriDiskGeometry::SigmaZ() const
{
    return 0.;
}

//////////////////////////////////////////////////////////////////////
