/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "ClipGeometryDecorator.hpp"
#include "FatalError.hpp"

////////////////////////////////////////////////////////////////////

void ClipGeometryDecorator::setupSelfAfter()
{
    Geometry::setupSelfAfter();

    // estimate the original geometry's mass in the removed region
    int Nsamples = 10000;
    int Ninside = 0;
    for (int k = 0; k < Nsamples; k++)
    {
        Position bfr = _geometry->generatePosition();
        if (inside(bfr)) Ninside++;
    }
    double chi = Ninside / (1.0 * Nsamples);
    if (_remove == Remove::Outside) chi = 1.0 - chi;
    if (chi > 0.9) throw FATALERROR("Clip decorator removes more than 90% of the original mass");
    _norm = 1.0 / (1.0 - chi);
}

//////////////////////////////////////////////////////////////////////

double ClipGeometryDecorator::density(Position bfr) const
{
    bool removed = inside(bfr);
    if (_remove == Remove::Outside) removed = !removed;
    return removed ? 0.0 : _geometry->density(bfr) * _norm;
}

////////////////////////////////////////////////////////////////////

Position ClipGeometryDecorator::generatePosition() const
{
    while (true)
    {
        Position bfr = _geometry->generatePosition();
        bool removed = inside(bfr);
        if (_remove == Remove::Outside) removed = !removed;
        if (!removed) return bfr;
    }
}

////////////////////////////////////////////////////////////////////

double ClipGeometryDecorator::SigmaX() const
{
    return _geometry->SigmaX();
}

////////////////////////////////////////////////////////////////////

double ClipGeometryDecorator::SigmaY() const
{
    return _geometry->SigmaY();
}

////////////////////////////////////////////////////////////////////

double ClipGeometryDecorator::SigmaZ() const
{
    return _geometry->SigmaZ();
}

////////////////////////////////////////////////////////////////////

double ClipGeometryDecorator::norm() const
{
    return _norm;
}

////////////////////////////////////////////////////////////////////
