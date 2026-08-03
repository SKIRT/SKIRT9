/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "SpheroidalGeometryDecorator.hpp"
#include "Random.hpp"

////////////////////////////////////////////////////////////////////

double SpheroidalGeometryDecorator::density(double R, double z) const
{
    double m = sqrt(R * R + z * z / (_q * _q));
    return 1.0 / _q * _geometry->density(m);
}

////////////////////////////////////////////////////////////////////

Position SpheroidalGeometryDecorator::generatePosition() const
{
    Position bfrs = _geometry->generatePosition();
    double xs, ys, zs;
    bfrs.cartesian(xs, ys, zs);
    return Position(xs, ys, _q * zs);
}

////////////////////////////////////////////////////////////////////

double SpheroidalGeometryDecorator::SigmaR() const
{
    return 1.0 / _q * _geometry->Sigmar();
}

////////////////////////////////////////////////////////////////////

double SpheroidalGeometryDecorator::SigmaZ() const
{
    return 2.0 * _geometry->Sigmar();
}

////////////////////////////////////////////////////////////////////
