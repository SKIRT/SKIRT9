/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "CombineGeometryDecorator.hpp"
#include "Random.hpp"

////////////////////////////////////////////////////////////////////

void CombineGeometryDecorator::setupSelfBefore()
{
    Geometry::setupSelfBefore();

    double sum = _firstWeight + _secondWeight;
    _w1 = _firstWeight / sum;
    _w2 = _secondWeight / sum;
}

//////////////////////////////////////////////////////////////////////

int CombineGeometryDecorator::dimension() const
{
    return max(_firstGeometry->dimension(), _secondGeometry->dimension());
}

////////////////////////////////////////////////////////////////////

double CombineGeometryDecorator::density(Position bfr) const
{
    return _w1 * _firstGeometry->density(bfr) + _w2 * _secondGeometry->density(bfr);
}

////////////////////////////////////////////////////////////////////

Position CombineGeometryDecorator::generatePosition() const
{
    double X = random()->uniform();
    if (X < _w1)
        return _firstGeometry->generatePosition();
    else
        return _secondGeometry->generatePosition();
}

////////////////////////////////////////////////////////////////////

double CombineGeometryDecorator::SigmaX() const
{
    return _w1 * _firstGeometry->SigmaX() + _w2 * _secondGeometry->SigmaX();
}

////////////////////////////////////////////////////////////////////

double CombineGeometryDecorator::SigmaY() const
{
    return _w1 * _firstGeometry->SigmaY() + _w2 * _secondGeometry->SigmaY();
}

////////////////////////////////////////////////////////////////////

double CombineGeometryDecorator::SigmaZ() const
{
    return _w1 * _firstGeometry->SigmaZ() + _w2 * _secondGeometry->SigmaZ();
}

////////////////////////////////////////////////////////////////////
