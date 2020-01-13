/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "HyperboloidShellGeometry.hpp"
#include "FatalError.hpp"
#include "Random.hpp"

//////////////////////////////////////////////////////////////////////

void HyperboloidShellGeometry::setupSelfBefore()
{
    AxGeometry::setupSelfBefore();

    // verify property values
    if (_DeltaOut <= _DeltaIn) throw FATALERROR("the outer opening angle should be larger than the inner angle");
    if (_aout <= _ain) throw FATALERROR("the outer real axis should be larger than the inner axis");

    // determine the maximal vertical extent of the hyperboloid
    _zmax = _Dout * cos(_DeltaOut);

    // determine the radius of the outer wall top cross section
    _bout = _Dout * sin(_DeltaOut);

    // determine the imaginary axis of the hyperboloid outer wall
    _cout = _aout * _zmax / sqrt(_bout * _bout - _aout * _aout);

    // determine the radial extent of the hyperboloid inner wall
    _Din = _zmax / cos(_DeltaIn);

    // determine the radius of the inner wall top cross section
    _bin = _Din * sin(_DeltaIn);

    // determine the imaginary axis of the hyperboloid inner wall
    _cin = _ain * _zmax / sqrt(_bin * _bin - _ain * _ain);

    // determine the normalization factor
    _A = 3. / 2. / M_PI / _zmax / (2 * _aout * _aout + _bout * _bout - 2 * _ain * _ain - _bin * _bin);
}

///////////////////////////////////////////////////////////////////

double HyperboloidShellGeometry::density(double R, double z) const
{
    if (z > _zmax) return 0.;
    if (z < -_zmax) return 0.;
    if (R < _ain) return 0.;
    if ((z >= 0.) && (z < _cout / _aout * sqrt(R * R - _aout * _aout))) return 0.;
    if ((z >= 0.) && (z > _cin / _ain * sqrt(R * R - _ain * _ain))) return 0.;
    if ((z < 0.) && (z > -_cout / _aout * sqrt(R * R - _aout * _aout))) return 0.;
    if ((z < 0.) && (z < -_cin / _ain * sqrt(R * R - _ain * _ain))) return 0.;
    return _A;
}

//////////////////////////////////////////////////////////////////////

Position HyperboloidShellGeometry::generatePosition() const
{
    while (true)
    {
        double R = _bout * sqrt(random()->uniform());
        double z = _zmax * (2. * random()->uniform() - 1.);
        if (density(R, z))
        {
            double phi = 2. * M_PI * random()->uniform();
            return Position(R, phi, z, Position::CoordinateSystem::CYLINDRICAL);
        }
    }
}

//////////////////////////////////////////////////////////////////////

double HyperboloidShellGeometry::SigmaR() const
{
    return _A * (_aout - _ain);
}

//////////////////////////////////////////////////////////////////////

double HyperboloidShellGeometry::SigmaZ() const
{
    return 0.;
}

//////////////////////////////////////////////////////////////////////
