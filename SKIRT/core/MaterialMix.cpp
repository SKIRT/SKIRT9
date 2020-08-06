/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "MaterialMix.hpp"
#include "FatalError.hpp"
#include "Random.hpp"

////////////////////////////////////////////////////////////////////

void MaterialMix::setupSelfBefore()
{
    SimulationItem::setupSelfBefore();

    _random = find<Random>();
}

////////////////////////////////////////////////////////////////////

bool MaterialMix::hasPolarizedScattering() const
{
    return false;
}

////////////////////////////////////////////////////////////////////

bool MaterialMix::hasPolarizedAbsorption() const
{
    return false;
}

////////////////////////////////////////////////////////////////////

bool MaterialMix::hasPolarizedEmission() const
{
    return false;
}

////////////////////////////////////////////////////////////////////

bool MaterialMix::hasResonantScattering() const
{
    return false;
}

////////////////////////////////////////////////////////////////////

bool MaterialMix::hasStochasticDustEmission() const
{
    return false;
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
    return 2. * random()->uniform() - 1.;
}

////////////////////////////////////////////////////////////////////

double MaterialMix::phaseFunctionValue(double /*lambda*/, double /*theta*/, double /*phi*/,
                                       const StokesVector* /*sv*/) const
{
    throw FATALERROR("This function implementation should never be called");
}

////////////////////////////////////////////////////////////////////

std::pair<double, double> MaterialMix::generateAnglesFromPhaseFunction(double /*lambda*/,
                                                                       const StokesVector* /*sv*/) const
{
    throw FATALERROR("This function implementation should never be called");
}

////////////////////////////////////////////////////////////////////

void MaterialMix::applyMueller(double /*lambda*/, double /*theta*/, StokesVector* /*sv*/) const
{
    throw FATALERROR("This function implementation should never be called");
}

////////////////////////////////////////////////////////////////////

const Array& MaterialMix::thetaGrid() const
{
    throw FATALERROR("This function implementation should never be called");
}

////////////////////////////////////////////////////////////////////

const Array& MaterialMix::sectionsAbs(double /*lambda*/) const
{
    throw FATALERROR("This function implementation should never be called");
}

////////////////////////////////////////////////////////////////////

const Array& MaterialMix::sectionsAbspol(double /*lambda*/) const
{
    throw FATALERROR("This function implementation should never be called");
}

////////////////////////////////////////////////////////////////////
