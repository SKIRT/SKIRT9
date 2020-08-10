/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "LyaNeutralHydrogenMaterialMix.hpp"
#include "Configuration.hpp"
#include "Constants.hpp"
#include "LyaUtils.hpp"
#include "MediumState.hpp"
#include "PhotonPacket.hpp"
#include "Random.hpp"
#include "StokesVector.hpp"

////////////////////////////////////////////////////////////////////

void LyaNeutralHydrogenMaterialMix::setupSelfBefore()
{
    MaterialMix::setupSelfBefore();

    _dpf.initialize(random(), includePolarization());
}

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

void LyaNeutralHydrogenMaterialMix::performScattering(const MediumState* state, PhotonPacket* pp) const
{
    // calculate the perceived wavelength in the cell
    Vec bfv = state->bulkVelocity();
    double lambda = config()->hasMovingMedia()
                        ? pp->perceivedWavelength(bfv, config()->lyaExpansionRate() * pp->interactionDistance())
                        : pp->wavelength();

    // draw a random atom velocity & phase function, unless a peel-off stored this already
    if (!pp->hasLyaScatteringInfo())
    {
        double T = state->temperature();
        double nH = state->numberDensity();
        pp->setLyaScatteringInfo(LyaUtils::sampleAtomVelocity(lambda, T, nH, pp->direction(), config(), random()));
    }

    // draw the outgoing direction from the dipole or the isotropic phase function
    // and, if required, update the polarization state of the photon packet
    Direction bfknew;
    if (pp->lyaDipole())
    {
        bfknew = _dpf.performScattering(pp->direction(), pp);
    }
    else
    {
        bfknew = random()->direction();
        if (includePolarization()) pp->setUnpolarized();
    }

    // Doppler-shift the photon packet wavelength into and out of the atom frame
    lambda = LyaUtils::shiftWavelength(lambda, pp->lyaAtomVelocity(), pp->direction(), bfknew);

    // if the material in the cell has a nonzero bulk velocity, determine the Doppler-shifted outgoing wavelength
    if (!bfv.isNull()) lambda = PhotonPacket::shiftedEmissionWavelength(lambda, bfknew, bfv);

    // set the scattering event in the photon packet
    pp->scatter(bfknew, lambda);
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
