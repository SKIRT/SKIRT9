/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "AxAngularDistribution.hpp"
#include "FatalError.hpp"
#include "Random.hpp"

//////////////////////////////////////////////////////////////////////

void AxAngularDistribution::setupSelfBefore()
{
    AngularDistribution::setupSelfBefore();

    Vec sym(symmetryX(), symmetryY(), symmetryZ());
    if (sym.norm() == 0) throw FATALERROR("Symmetry axis direction cannot be null vector");
    _sym = Direction(sym / sym.norm());
}

//////////////////////////////////////////////////////////////////////

int AxAngularDistribution::dimension() const
{
    return symmetryX() || symmetryY() ? 3 : 2;
}

//////////////////////////////////////////////////////////////////////

double AxAngularDistribution::probabilityForDirection(Direction bfk) const
{
    return probabilityForInclinationCosine(Vec::dot(_sym, bfk));
}

//////////////////////////////////////////////////////////////////////

Direction AxAngularDistribution::generateDirection() const
{
    return random()->direction(_sym, generateInclinationCosine());
}

//////////////////////////////////////////////////////////////////////
