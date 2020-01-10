/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "NetzerAngularDistribution.hpp"
#include "NR.hpp"
#include "Random.hpp"

//////////////////////////////////////////////////////////////////////

void NetzerAngularDistribution::setupSelfBefore()
{
    AxAngularDistribution::setupSelfBefore();

    // grid with values of the cumulative luminosity distribution
    int n = 400;
    NR::buildLinearGrid(_costhetav, -1., +1., n);
    _Xv.resize(n + 1);
    _Xv[0] = 0;
    for (int i = 1; i < n; i++)
    {
        double ct = _costhetav[i];
        double sign = ct > 0 ? 1. : -1;
        _Xv[i] = (1. / 2.) + (2. / 7.) * ct * ct * ct + sign * (3. / 14.) * ct * ct;
    }
    _Xv[n] = 1.;
}

//////////////////////////////////////////////////////////////////////

double NetzerAngularDistribution::probabilityForInclinationCosine(double costheta) const
{
    double sign = costheta > 0 ? 1. : -1;
    return (6. / 7.) * costheta * (2. * costheta + sign);
}

//////////////////////////////////////////////////////////////////////

double NetzerAngularDistribution::generateInclinationCosine() const
{
    return random()->cdfLinLin(_costhetav, _Xv);
}

//////////////////////////////////////////////////////////////////////
