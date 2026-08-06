/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "ConicalShellGeometry.hpp"
#include "FatalError.hpp"
#include "Random.hpp"
#include "SpecialFunctions.hpp"

//////////////////////////////////////////////////////////////////////

void ConicalShellGeometry::setupSelfBefore()
{
    AxGeometry::setupSelfBefore();

    // verify property values
    if (_DeltaOut <= _DeltaIn) throw FATALERROR("the outer angle of the shell should be larger than the inner angle");
    if (_rmax <= _rmin) throw FATALERROR("The maximum radius of the shell should be larger than the minimum radius");

    // cache frequently used values
    _sinDeltaIn = sin(_DeltaIn);
    _sinDeltaOut = sin(_DeltaOut);
    _cosDelta = cos((_DeltaOut + _DeltaIn) / 2.);
    _smin = SpecialFunctions::gln(_p - 2.0, _rmin);
    _sdiff = SpecialFunctions::gln2(_p - 2.0, _rmax, _rmin);
    _tmin = pow(_rmin, 3.0 - _p);
    _tmax = pow(_rmax, 3.0 - _p);

    // determine the normalization factor
    if (_q > 1e-3)
        _A = _q * 0.25 / M_PI / _sdiff / (exp(-_q * _sinDeltaIn) - exp(-_q * _sinDeltaOut));
    else
        _A = 0.25 / M_PI / _sdiff / (_sinDeltaOut - _sinDeltaIn);
}

//////////////////////////////////////////////////////////////////////

double ConicalShellGeometry::density(double R, double z) const
{
    double r = sqrt(R * R + z * z);
    double costheta = z / r;

    if (r >= _rmax) return 0.0;
    if (_rani)
    {
        double rminani = _rmin * sqrt(6. / 7. * fabs(costheta) * (2. * fabs(costheta) + 1));
        if (r <= rminani || r < _rcut) return 0.0;
    }
    else
    {
        if (r <= _rmin) return 0.0;
    }

    if (fabs(costheta) >= _sinDeltaOut) return 0.0;
    if (fabs(costheta) <= _sinDeltaIn) return 0.0;
    return _A * pow(r, -_p) * exp(-_q * fabs(costheta));
}

//////////////////////////////////////////////////////////////////////

Position ConicalShellGeometry::generatePosition() const
{
    while (true)
    {
        double X = random()->uniform();
        double r = 0.0;
        if (fabs(_p - 3.0) < 1e-2)
        {
            double s = _smin + X * _sdiff;
            r = SpecialFunctions::gexp(_p - 2.0, s);
        }
        else
        {
            double z = (1.0 - X) * _tmin + X * _tmax;
            r = pow(z, 1.0 / (3.0 - _p));
        }

        X = random()->uniform();
        double costheta = 0.0;
        if (_q < 1e-3)
            costheta = (1.0 - 2.0 * X) * _sinDeltaOut;
        else
        {
            double B = 1.0 - exp(-_q * _sinDeltaOut);
            costheta = (X < 0.5) ? -log(1.0 - B * (1.0 - 2.0 * X)) / _q : log(1.0 - B * (2.0 * X - 1.0)) / _q;
        }
        double theta = acos(costheta);
        X = random()->uniform();
        double phi = 2.0 * M_PI * X;
        Position bfr(r, theta, phi, Position::CoordinateSystem::SPHERICAL);
        if (density(bfr.cylRadius(), bfr.height())) return bfr;
    }
}

//////////////////////////////////////////////////////////////////////

double ConicalShellGeometry::SigmaR() const
{
    return _A * exp(-_q * _cosDelta) * SpecialFunctions::gln2(_p, _rmax, _rmin);
}

//////////////////////////////////////////////////////////////////////

double ConicalShellGeometry::SigmaZ() const
{
    return 0.0;
}

//////////////////////////////////////////////////////////////////////
