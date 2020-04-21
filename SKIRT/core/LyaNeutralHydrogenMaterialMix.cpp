/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "LyaNeutralHydrogenMaterialMix.hpp"
#include "Constants.hpp"
#include "LyaUtils.hpp"
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
