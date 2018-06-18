/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "NetzerAngularDistribution.hpp"
#include "Random.hpp"

//////////////////////////////////////////////////////////////////////

void NetzerAngularDistribution::setupSelfBefore()
{
    AngularDistribution::setupSelfBefore();

    // grid with values of the cumulative luminosity distribution
    int N = 401;
    _thetav.resize(N);
    _Xv.resize(N);
    _thetav[0] = 0;
    _Xv[0] = 0;
    for (int i=1; i<N-1; i++)
    {
        _thetav[i] = M_PI*i/(N-1);
        double ct = cos(_thetav[i]);
        double sign = ct > 0 ? 1. : -1;
        _Xv[i] = (1./2.) - (2./7.)*ct*ct*ct - sign*(3./14.)*ct*ct;
    }
    _thetav[N-1] = M_PI;
    _Xv[N-1] = 1.0;
}

//////////////////////////////////////////////////////////////////////

int NetzerAngularDistribution::dimension() const
{
    return 2;
}

//////////////////////////////////////////////////////////////////////

double NetzerAngularDistribution::probabilityForDirection(Direction bfk) const
{
    double theta, phi;
    bfk.spherical(theta, phi);
    double ct = cos(theta);
    double sign = ct > 0 ? 1. : -1;
    return (6./7.) * ct * (2.*ct + sign);
}

//////////////////////////////////////////////////////////////////////

Direction NetzerAngularDistribution::generateDirection() const
{
    double theta = random()->cdfLinLin(_thetav,_Xv);
    double phi = 2.0*M_PI*random()->uniform();
    return Direction(theta,phi);
}

//////////////////////////////////////////////////////////////////////
