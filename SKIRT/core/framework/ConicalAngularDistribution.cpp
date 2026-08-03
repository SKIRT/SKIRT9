/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "ConicalAngularDistribution.hpp"
#include "Random.hpp"

////////////////////////////////////////////////////////////////////

void ConicalAngularDistribution::setupSelfBefore()
{
    AxAngularDistribution::setupSelfBefore();
    _cosDelta = cos(_openingAngle);
}

////////////////////////////////////////////////////////////////////

double ConicalAngularDistribution::probabilityForInclinationCosine(double costheta) const
{
    if (abs(costheta) > _cosDelta)
        return 1.0 / (1.0 - _cosDelta);
    else
        return 0.;
}

////////////////////////////////////////////////////////////////////

double ConicalAngularDistribution::generateInclinationCosine() const
{
    double X = random()->uniform();
    if (X < 0.5)
        return 1.0 - 2.0 * X * (1.0 - _cosDelta);
    else
        return 1.0 - 2.0 * _cosDelta - 2.0 * X * (1.0 - _cosDelta);
}

////////////////////////////////////////////////////////////////////
