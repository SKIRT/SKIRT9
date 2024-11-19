
/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "RedistributeGeometryDecorator.hpp"
#include "Random.hpp"

////////////////////////////////////////////////////////////////////

void RedistributeGeometryDecorator::setupSelfAfter()
{
    Geometry::setupSelfAfter();

    int Nsamples = 10000;
    double sum = 0.;
    for (int k = 0; k < Nsamples; k++)
    {
        Position bfr = _geometry->generatePosition();
        sum += weight(bfr);
    }
    _norm = Nsamples / sum;
    _maxWeight = maxWeight();
}

//////////////////////////////////////////////////////////////////////

double RedistributeGeometryDecorator::density(Position bfr) const
{
    return _norm * _geometry->density(bfr) * weight(bfr);
}

////////////////////////////////////////////////////////////////////

Position RedistributeGeometryDecorator::generatePosition() const
{
    while (true)
    {
        Position bfr = _geometry->generatePosition();
        double t = random()->uniform() * _maxWeight / weight(bfr);
        if (t <= 1.) return bfr;
    }
}

////////////////////////////////////////////////////////////////////

double RedistributeGeometryDecorator::SigmaX() const
{
    return _geometry->SigmaX();
}

////////////////////////////////////////////////////////////////////

double RedistributeGeometryDecorator::SigmaY() const
{
    return _geometry->SigmaY();
}

////////////////////////////////////////////////////////////////////

double RedistributeGeometryDecorator::SigmaZ() const
{
    return _geometry->SigmaZ();
}

////////////////////////////////////////////////////////////////////
