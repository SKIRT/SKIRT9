/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "LyaNeutralHydrogenMaterialMix.hpp"
#include "Constants.hpp"
#include "LyaUtils.hpp"
#include "MediumState.hpp"
#include "PhotonPacket.hpp"
#include "Random.hpp"
#include "StokesVector.hpp"

////////////////////////////////////////////////////////////////////

MaterialMix::MaterialType LyaNeutralHydrogenMaterialMix::materialType() const
{
    return MaterialType::Gas;
}

////////////////////////////////////////////////////////////////////

MaterialMix::ScatteringMode LyaNeutralHydrogenMaterialMix::scatteringMode() const
{
    return includePolarization() ? ScatteringMode::LyaPolarization : ScatteringMode::Lya;
}

////////////////////////////////////////////////////////////////////

bool LyaNeutralHydrogenMaterialMix::hasPolarizedScattering() const
{
    return includePolarization();
}

////////////////////////////////////////////////////////////////////

bool LyaNeutralHydrogenMaterialMix::hasResonantScattering() const
{
    return true;
}

////////////////////////////////////////////////////////////////////

double LyaNeutralHydrogenMaterialMix::mass() const
{
    return Constants::Mproton();
}

////////////////////////////////////////////////////////////////////

double LyaNeutralHydrogenMaterialMix::sectionAbs(double /*lambda*/) const
{
    return 0.;
}

////////////////////////////////////////////////////////////////////

double LyaNeutralHydrogenMaterialMix::sectionSca(double lambda) const
{
    return LyaUtils::section(lambda, defaultTemperature());
}

////////////////////////////////////////////////////////////////////

double LyaNeutralHydrogenMaterialMix::sectionExt(double lambda) const
{
    return LyaUtils::section(lambda, defaultTemperature());
}

////////////////////////////////////////////////////////////////////

double LyaNeutralHydrogenMaterialMix::opacityAbs(double /*lambda*/, const MediumState* /*state*/,
                                                 const PhotonPacket* /*pp*/) const
{
    return 0.;
}

////////////////////////////////////////////////////////////////////

double LyaNeutralHydrogenMaterialMix::opacitySca(double lambda, const MediumState* state,
                                                 const PhotonPacket* /*pp*/) const
{
    double n = state->numberDensity();
    return n > 0. ? n * LyaUtils::section(lambda, state->temperature()) : 0.;
}

////////////////////////////////////////////////////////////////////

double LyaNeutralHydrogenMaterialMix::opacityExt(double lambda, const MediumState* state,
                                                 const PhotonPacket* /*pp*/) const
{
    double n = state->numberDensity();
    return n > 0. ? n * LyaUtils::section(lambda, state->temperature()) : 0.;
}

////////////////////////////////////////////////////////////////////

double LyaNeutralHydrogenMaterialMix::equilibriumTemperature(const Array& /*Jv*/) const
{
    return defaultTemperature();
}

////////////////////////////////////////////////////////////////////

Array LyaNeutralHydrogenMaterialMix::emissivity(const Array& /*Jv*/) const
{
    return Array();
}

////////////////////////////////////////////////////////////////////
