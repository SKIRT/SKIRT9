/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "Direction.hpp"
#include "FatalError.hpp"

//////////////////////////////////////////////////////////////////////

Direction::Direction(double theta, double phi)
{
    static const double eps = 1e-8;
    if (theta < -eps || theta > M_PI + eps)
    {
        throw FATALERROR("Theta should be between 0 and pi.");
    }
    else if (theta <= eps)
    {
        _x = 0;
        _y = 0;
        _z = 1;
    }
    else if (theta >= M_PI - eps)
    {
        _x = 0;
        _y = 0;
        _z = -1;
    }
    else
    {
        double sintheta = sin(theta);
        _x = sintheta * cos(phi);
        _y = sintheta * sin(phi);
        _z = cos(theta);
    }
}

//////////////////////////////////////////////////////////////////////

void Direction::cartesian(double& kx, double& ky, double& kz) const
{
    kx = _x;
    ky = _y;
    kz = _z;
}

//////////////////////////////////////////////////////////////////////

void Direction::spherical(double& theta, double& phi) const
{
    theta = acos(_z);
    if (_x == 0 && _y == 0)
        phi = 0;
    else
        phi = atan2(_y, _x);
}

//////////////////////////////////////////////////////////////////////
