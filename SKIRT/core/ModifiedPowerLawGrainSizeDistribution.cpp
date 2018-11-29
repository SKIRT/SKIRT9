/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "ModifiedPowerLawGrainSizeDistribution.hpp"

////////////////////////////////////////////////////////////////////

double ModifiedPowerLawGrainSizeDistribution::dnda(double a) const
{
    return    pow(a,_alpha)
            * pow(1. + fabs(_zeta)*pow(a/_au,_eta), _zeta>=0 ? 1 : -1)  // curvature
            * ( a<=_at ? 1. : exp(-pow((a-_at)/_ac,_gamma)) );          // exponential decay
}

////////////////////////////////////////////////////////////////////
