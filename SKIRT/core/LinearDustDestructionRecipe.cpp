/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "LinearDustDestructionRecipe.hpp"
#include "FatalError.hpp"

////////////////////////////////////////////////////////////////////

void LinearDustDestructionRecipe::setupSelfBefore()
{
    DustDestructionRecipe::setupSelfBefore();

    if (maxGraphiteTemperature() < minGraphiteTemperature() || maxSilicateTemperature() < minSilicateTemperature())
        throw FATALERROR("Maximum temperature is below minimum temperature");
}

////////////////////////////////////////////////////////////////////

double LinearDustDestructionRecipe::densityFraction(bool graphite, double /*a*/, const Array& /*Jv*/, double T) const
{
    // handle the cases outside the configured temperature interval
    double Tmin = graphite ? minGraphiteTemperature() : minSilicateTemperature();
    if (T <= Tmin) return 1.;
    double Tmax = graphite ? maxGraphiteTemperature() : maxSilicateTemperature();
    if (T >= Tmax) return 0.;

    // we now know that Tmin < T < Tmax and thus that Tmax - Tmin > 0
    return (Tmax - T) / (Tmax - Tmin);
}

////////////////////////////////////////////////////////////////////
