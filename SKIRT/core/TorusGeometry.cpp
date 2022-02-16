/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "TorusGeometry.hpp"
#include "Random.hpp"
#include "SpecialFunctions.hpp"

//////////////////////////////////////////////////////////////////////

void TorusGeometry::setupSelfBefore()
{
    AxGeometry::setupSelfBefore();

    // cache frequently used values
    _sinDelta = sin(_Delta);
    _smin = SpecialFunctions::gln(_p - 2.0, _rmin);
    _sdiff = SpecialFunctions::gln2(_p - 2.0, _rmax, _rmin);
    _tmin = pow(_rmin, 3.0 - _p);
    _tmax = pow(_rmax, 3.0 - _p);

    // determine the normalization factor
    if (_q > 1e-3)
        _A = _q * 0.25 / M_PI / _sdiff / (1.0 - exp(-_q * _sinDelta));
    else
        _A = 0.25 / M_PI / _sdiff / _sinDelta;
}

//////////////////////////////////////////////////////////////////////

double TorusGeometry::density(double R, double z) const
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

    if (fabs(costheta) >= _sinDelta) return 0.0;
    return _A * pow(r, -_p) * exp(-_q * fabs(costheta));
}

//////////////////////////////////////////////////////////////////////

Position TorusGeometry::generatePosition() const
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
            costheta = (1.0 - 2.0 * X) * _sinDelta;
        else
        {
            double B = 1.0 - exp(-_q * _sinDelta);
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

double TorusGeometry::SigmaR() const
{
    return _A * SpecialFunctions::gln2(_p, _rmax, _rmin);
}

//////////////////////////////////////////////////////////////////////

double TorusGeometry::SigmaZ() const
{
    return 0.0;
}

//////////////////////////////////////////////////////////////////////
