/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "SineSquarePolarizationProfile.hpp"
#include "FatalError.hpp"

////////////////////////////////////////////////////////////////////

void SineSquarePolarizationProfile::setupSelfBefore()
{
    PolarizationProfile::setupSelfBefore();

    Vec sym(symmetryX(), symmetryY(), symmetryZ());
    if (sym.norm() == 0) throw FATALERROR("Symmetry axis direction cannot be null vector");
    _sym = Direction(sym / sym.norm());

    _cos2gamma = cos(2 * polarizationAngle());
    _sin2gamma = sin(2 * polarizationAngle());
}

////////////////////////////////////////////////////////////////////

int SineSquarePolarizationProfile::dimension() const
{
    return symmetryX() || symmetryY() ? 3 : 2;
}

////////////////////////////////////////////////////////////////////

StokesVector SineSquarePolarizationProfile::polarizationForDirection(Direction bfk) const
{
    double costheta = Vec::dot(_sym, bfk);
    if (abs(costheta) > 0.99999)
    {
        return StokesVector();
    }
    else
    {
        double PL = maxPolarizationDegree() * (1 - costheta) * (1 + costheta);
        double Q = PL * _cos2gamma;
        double U = PL * _sin2gamma;
        Direction n(Vec::cross(_sym, bfk));
        return StokesVector(1., Q, U, 0., n);
    }
}

////////////////////////////////////////////////////////////////////
