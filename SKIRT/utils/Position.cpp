/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "Position.hpp"

//////////////////////////////////////////////////////////////////////

Position::Position(double u, double v, double w, CoordinateSystem coordtype) : Vec()
{
    switch (coordtype)
    {
        case CoordinateSystem::CARTESIAN:
        {
            _x = u;
            _y = v;
            _z = w;
            break;
        }
        case CoordinateSystem::CYLINDRICAL:
        {
            double R = u;
            double phi = v;
            _x = R * cos(phi);
            _y = R * sin(phi);
            _z = w;
            break;
        }
        case CoordinateSystem::SPHERICAL:
        {
            double r = u;
            double theta = v;
            double phi = w;
            double costheta = cos(theta);
            double sintheta = sin(theta);
            double cosphi = cos(phi);
            double sinphi = sin(phi);
            _x = r * sintheta * cosphi;
            _y = r * sintheta * sinphi;
            _z = r * costheta;
            break;
        }
    }
}

//////////////////////////////////////////////////////////////////////

Position::Position(double r, Direction bfk) : Vec(r * bfk) {}

//////////////////////////////////////////////////////////////////////

Position::Position(Direction bfk) : Vec(bfk) {}

//////////////////////////////////////////////////////////////////////

double Position::radius() const
{
    return norm();
}

//////////////////////////////////////////////////////////////////////

double Position::cylRadius() const
{
    return sqrt(_x * _x + _y * _y);
}

//////////////////////////////////////////////////////////////////////

double Position::height() const
{
    return _z;
}

//////////////////////////////////////////////////////////////////////

void Position::cartesian(double& x, double& y, double& z) const
{
    x = _x;
    y = _y;
    z = _z;
}

//////////////////////////////////////////////////////////////////////

void Position::spherical(double& r, double& theta, double& phi) const
{
    r = radius();
    if (r == 0)
    {
        theta = 0;
        phi = 0;
    }
    else
    {
        theta = acos(_z / r);
        phi = atan2(_y, _x);
    }
}

//////////////////////////////////////////////////////////////////////

void Position::cylindrical(double& R, double& phi, double& z) const
{
    R = cylRadius();
    phi = atan2(_y, _x);
    z = _z;
}

//////////////////////////////////////////////////////////////////////
