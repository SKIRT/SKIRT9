/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "OffsetGeometryDecorator.hpp"

////////////////////////////////////////////////////////////////////

int OffsetGeometryDecorator::dimension() const
{
    return (_offsetX || _offsetY || _geometry->dimension() == 3) ? 3 : 2;
}

////////////////////////////////////////////////////////////////////

double OffsetGeometryDecorator::density(Position bfr) const
{
    double x, y, z;
    bfr.cartesian(x, y, z);
    return _geometry->density(Position(x - _offsetX, y - _offsetY, z - _offsetZ));
}

////////////////////////////////////////////////////////////////////

Position OffsetGeometryDecorator::generatePosition() const
{
    Position bfr = _geometry->generatePosition();
    double x, y, z;
    bfr.cartesian(x, y, z);
    return Position(x + _offsetX, y + _offsetY, z + _offsetZ);
}

////////////////////////////////////////////////////////////////////

double OffsetGeometryDecorator::SigmaX() const
{
    return _geometry->SigmaX();
}

////////////////////////////////////////////////////////////////////

double OffsetGeometryDecorator::SigmaY() const
{
    return _geometry->SigmaY();
}

////////////////////////////////////////////////////////////////////

double OffsetGeometryDecorator::SigmaZ() const
{
    return _geometry->SigmaZ();
}

////////////////////////////////////////////////////////////////////
