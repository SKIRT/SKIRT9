/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "ConfigurableDustMix.hpp"

//////////////////////////////////////////////////////////////////////

void ConfigurableDustMix::setupSelfBefore()
{
    MultiGrainDustMix::setupSelfBefore();

    for (auto population : _populations) addPopulation(population);
}

//////////////////////////////////////////////////////////////////////

DustMix::ScatteringMode ConfigurableDustMix::scatteringMode() const
{
    switch (scatteringType())
    {
        case ScatteringType::HenyeyGreenstein: return ScatteringMode::HenyeyGreenstein;
        case ScatteringType::MaterialPhaseFunction: return ScatteringMode::MaterialPhaseFunction;
        case ScatteringType::SphericalPolarization: return ScatteringMode::SphericalPolarization;
        case ScatteringType::SpheroidalPolarization: return ScatteringMode::SpheroidalPolarization;
    }
    return ScatteringMode::HenyeyGreenstein;  // to satisfy gcc compiler
}

//////////////////////////////////////////////////////////////////////

bool ConfigurableDustMix::hasPolarizedScattering() const
{
    return scatteringType() == ScatteringType::SphericalPolarization
           || scatteringType() == ScatteringType::SpheroidalPolarization;
}

//////////////////////////////////////////////////////////////////////

bool ConfigurableDustMix::hasPolarizedEmission() const
{
    return scatteringType() == ScatteringType::SpheroidalPolarization;
}

//////////////////////////////////////////////////////////////////////
