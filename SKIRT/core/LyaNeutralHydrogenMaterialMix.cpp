/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "LyaNeutralHydrogenMaterialMix.hpp"
#include "Configuration.hpp"
#include "Constants.hpp"
#include "LyaUtils.hpp"
#include "MaterialState.hpp"
#include "PhotonPacket.hpp"
#include "Random.hpp"

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

vector<StateVariable> LyaNeutralHydrogenMaterialMix::specificStateVariableInfo() const
{
    return vector<StateVariable>{StateVariable::numberDensity(), StateVariable::temperature()};
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

double LyaNeutralHydrogenMaterialMix::opacityAbs(double /*lambda*/, const MaterialState* /*state*/,
                                                 const PhotonPacket* /*pp*/) const
{
    return 0.;
}

////////////////////////////////////////////////////////////////////

double LyaNeutralHydrogenMaterialMix::opacitySca(double lambda, const MaterialState* state,
                                                 const PhotonPacket* /*pp*/) const
{
    double n = state->numberDensity();
    return n > 0. ? n * LyaUtils::section(lambda, state->temperature()) : 0.;
}

////////////////////////////////////////////////////////////////////

double LyaNeutralHydrogenMaterialMix::opacityExt(double lambda, const MaterialState* state,
                                                 const PhotonPacket* /*pp*/) const
{
    double n = state->numberDensity();
    return n > 0. ? n * LyaUtils::section(lambda, state->temperature()) : 0.;
}

////////////////////////////////////////////////////////////////////

double LyaNeutralHydrogenMaterialMix::indicativeTemperature(const Array& /*Jv*/) const
{
    return defaultTemperature();
}

////////////////////////////////////////////////////////////////////

void LyaNeutralHydrogenMaterialMix::peeloffScattering(double& I, double& Q, double& U, double& V, double& lambda,
                                                      double w, Direction bfkobs, Direction bfky,
                                                      const MaterialState* state, PhotonPacket* pp) const
{
    // draw a random atom velocity & phase function, unless a previous peel-off stored this already
    if (!pp->hasLyaScatteringInfo())
    {
        double T = state->temperature();
        double nH = state->numberDensity();
        pp->setLyaScatteringInfo(LyaUtils::sampleAtomVelocity(lambda, T, nH, pp->direction(), config(), random()));
    }

    // add the contribution to the Stokes vector components depending on scattering type
    if (pp->lyaDipole())
    {
        // contribution of dipole scattering with or without polarization
        _dpf.peeloffScattering(I, Q, U, V, w, pp->direction(), bfkobs, bfky, pp);
    }
    else
    {
        // isotropic scattering removes polarization,
        // so the contribution is trivially 1 (multiplied by the weight for this component)
        I += w;
    }

    // for a random fraction of the events governed by the relative Lya contribution,
    // Doppler-shift the photon packet wavelength into and out of the atom frame
    if (random()->uniform() <= w)
        lambda = LyaUtils::shiftWavelength(lambda, pp->lyaAtomVelocity(), pp->direction(), bfkobs);
}

////////////////////////////////////////////////////////////////////

void LyaNeutralHydrogenMaterialMix::performScattering(double lambda, const MaterialState* state, PhotonPacket* pp) const
{
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

    // execute the scattering event in the photon packet
    pp->scatter(bfknew, state->bulkVelocity(), lambda);
}

////////////////////////////////////////////////////////////////////
