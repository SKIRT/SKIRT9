/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "MaterialMix.hpp"
#include "Random.hpp"

////////////////////////////////////////////////////////////////////

void MaterialMix::setupSelfBefore()
{
    SimulationItem::setupSelfBefore();

    _random = find<Random>();
}

////////////////////////////////////////////////////////////////////

MaterialMix::ScatteringMode MaterialMix::scatteringMode() const
{
    return ScatteringMode::HenyeyGreenstein;
}

////////////////////////////////////////////////////////////////////

double MaterialMix::asymmpar(double /*lambda*/) const
{
    return 0.;
}

////////////////////////////////////////////////////////////////////

double MaterialMix::phaseFunctionValueForCosine(double /*lambda*/, double /*costheta*/) const
{
    return 1.;
}

////////////////////////////////////////////////////////////////////

double MaterialMix::generateCosineFromPhaseFunction(double /*lambda*/) const
{
    return 2.*random()->uniform() - 1.;
}

////////////////////////////////////////////////////////////////////

double MaterialMix::phaseFunctionValue(double /*lambda*/, double /*theta*/, double /*phi*/,
                                       const StokesVector* /*sv*/) const
{
    return 1.;
}

////////////////////////////////////////////////////////////////////

std::pair<double, double> MaterialMix::generateAnglesFromPhaseFunction(double /*lambda*/,
                                                                       const StokesVector* /*sv*/) const
{
    return std::make_pair(0.,0.);
}

////////////////////////////////////////////////////////////////////

void MaterialMix::applyMueller(double /*lambda*/, double /*theta*/, StokesVector* /*sv*/) const
{
}

////////////////////////////////////////////////////////////////////
