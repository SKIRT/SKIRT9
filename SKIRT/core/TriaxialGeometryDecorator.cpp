/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "TriaxialGeometryDecorator.hpp"
#include "Random.hpp"

////////////////////////////////////////////////////////////////////

double TriaxialGeometryDecorator::density(Position bfr) const
{
    double x, y, z;
    bfr.cartesian(x, y, z);
    Position bfrs(x, y / _p, z / _q);
    return 1.0 / _p / _q * _geometry->density(bfrs);
}

////////////////////////////////////////////////////////////////////

Position TriaxialGeometryDecorator::generatePosition() const
{
    Position bfrs = _geometry->generatePosition();
    double xs, ys, zs;
    bfrs.cartesian(xs, ys, zs);
    return Position(xs, _p * ys, _q * zs);
}

////////////////////////////////////////////////////////////////////

double TriaxialGeometryDecorator::SigmaX() const
{
    return 2.0 / (_p * _q) * _geometry->Sigmar();
}

////////////////////////////////////////////////////////////////////

double TriaxialGeometryDecorator::SigmaY() const
{
    return 2.0 / _q * _geometry->Sigmar();
}

////////////////////////////////////////////////////////////////////

double TriaxialGeometryDecorator::SigmaZ() const
{
    return 2.0 / _p * _geometry->Sigmar();
}

////////////////////////////////////////////////////////////////////
