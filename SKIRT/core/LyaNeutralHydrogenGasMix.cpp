/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "LyaNeutralHydrogenGasMix.hpp"
#include "Constants.hpp"
#include "LyUtils.hpp"
#include "MaterialState.hpp"
#include "PhotonPacket.hpp"
#include "Random.hpp"

////////////////////////////////////////////////////////////////////

namespace
{
    // the combined Lya1 and Lya2 for hydrogen
    constexpr double lyaA = Constants::EinsteinALya();
    constexpr double lya = Constants::lambdaLya();
    constexpr double g = 3.;
    constexpr double kB = Constants::k();
    constexpr double mp = Constants::Mproton();
}

////////////////////////////////////////////////////////////////////

void LyaNeutralHydrogenGasMix::setupSelfBefore()
{
    MaterialMix::setupSelfBefore();

    _dpf.initialize(random(), includePolarization());
}

////////////////////////////////////////////////////////////////////

MaterialMix::MaterialType LyaNeutralHydrogenGasMix::materialType() const
{
    return MaterialType::Gas;
}

////////////////////////////////////////////////////////////////////

bool LyaNeutralHydrogenGasMix::hasPolarizedScattering() const
{
    return includePolarization();
}

////////////////////////////////////////////////////////////////////

bool LyaNeutralHydrogenGasMix::hasResonantScattering() const
{
    return true;
}

////////////////////////////////////////////////////////////////////

bool LyaNeutralHydrogenGasMix::hasExtraSpecificState() const
{
    return true;
}

////////////////////////////////////////////////////////////////////

bool LyaNeutralHydrogenGasMix::hasScatteringDispersion() const
{
    return true;
}

////////////////////////////////////////////////////////////////////

vector<StateVariable> LyaNeutralHydrogenGasMix::specificStateVariableInfo() const
{
    return vector<StateVariable>{StateVariable::numberDensity(), StateVariable::temperature()};
}

////////////////////////////////////////////////////////////////////

void LyaNeutralHydrogenGasMix::initializeSpecificState(MaterialState* state, double /*metallicity*/, double temperature,
                                                       const Array& /*params*/) const
{
    // leave the temperature at zero if the cell does not contain any material for this component
    if (state->numberDensity() > 0.)
    {
        // if no temperature was imported, use default value
        if (temperature < 0) temperature = defaultTemperature();

        // make sure the temperature is at least the local universe CMB temperature
        state->setTemperature(max(Constants::Tcmb(), temperature));
    }
}

////////////////////////////////////////////////////////////////////

double LyaNeutralHydrogenGasMix::mass() const
{
    return Constants::Mproton();
}

////////////////////////////////////////////////////////////////////

double LyaNeutralHydrogenGasMix::section(double lambda, double T) const
{
    double vth = sqrt(2. * kB / mp * T);
    double a = lyaA * lya / 4. / M_PI / vth;
    return LyUtils::section(lambda, lya, vth, lyaA, a, g);
}

////////////////////////////////////////////////////////////////////

double LyaNeutralHydrogenGasMix::sectionAbs(double /*lambda*/) const
{
    return 0.;
}

////////////////////////////////////////////////////////////////////

double LyaNeutralHydrogenGasMix::sectionSca(double lambda) const
{
    return section(lambda, defaultTemperature());
}

////////////////////////////////////////////////////////////////////

double LyaNeutralHydrogenGasMix::sectionExt(double lambda) const
{
    return section(lambda, defaultTemperature());
}

////////////////////////////////////////////////////////////////////

double LyaNeutralHydrogenGasMix::opacityAbs(double /*lambda*/, const MaterialState* /*state*/,
                                            const PhotonPacket* /*pp*/) const
{
    return 0.;
}

////////////////////////////////////////////////////////////////////

double LyaNeutralHydrogenGasMix::opacitySca(double lambda, const MaterialState* state, const PhotonPacket* /*pp*/) const
{
    double n = state->numberDensity();
    return n > 0. ? n * section(lambda, state->temperature()) : 0.;
}

////////////////////////////////////////////////////////////////////

double LyaNeutralHydrogenGasMix::opacityExt(double lambda, const MaterialState* state, const PhotonPacket* /*pp*/) const
{
    double n = state->numberDensity();
    return n > 0. ? n * section(lambda, state->temperature()) : 0.;
}

////////////////////////////////////////////////////////////////////

void LyaNeutralHydrogenGasMix::setScatteringInfoIfNeeded(PhotonPacket* pp, const MaterialState* state,
                                                         const double lambda) const
{
    auto scatinfo = pp->getScatteringInfo();
    if (!scatinfo->valid)
    {
        scatinfo->valid = true;

        double T = state->temperature();
        double nH = state->numberDensity();

        double vth = sqrt(2. * kB / mp * T);
        double a = lyaA * lya / 4. / M_PI / vth;
        double x = (lya - lambda) / lambda * Constants::c() / vth;

        // select the isotropic or the dipole phase function:
        // all wing events and 1/3 of core events are dipole, and the remaining 2/3 core events are isotropic,
        // where x=0.2 (in the atom frame) defines the transition between core and wings
        scatinfo->dipole = abs(x) > 0.2 || random()->uniform() < 1. / 3.;

        scatinfo->velocity =
            LyUtils::sampleAtomVelocity(lambda, lya, vth, a, T, nH, pp->direction(), config(), random());
    }
}

////////////////////////////////////////////////////////////////////

bool LyaNeutralHydrogenGasMix::peeloffScattering(double& I, double& Q, double& U, double& V, double& lambda,
                                                 Direction bfkobs, Direction bfky, const MaterialState* state,
                                                 const PhotonPacket* pp) const
{
    setScatteringInfoIfNeeded(const_cast<PhotonPacket*>(pp), state, lambda);
    auto scatinfo = const_cast<PhotonPacket*>(pp)->getScatteringInfo();

    // add the contribution to the Stokes vector components depending on scattering type
    if (scatinfo->dipole)
    {
        // contribution of dipole scattering with or without polarization
        _dpf.peeloffScattering(I, Q, U, V, pp->direction(), bfkobs, bfky, pp);
    }
    else
    {
        // isotropic scattering removes polarization, so the contribution is trivially 1
        I += 1.;
    }

    // Doppler-shift the photon packet wavelength into and out of the atom frame
    lambda = LyUtils::shiftWavelength(lambda, scatinfo->velocity, pp->direction(), bfkobs);

    return false;
}

////////////////////////////////////////////////////////////////////

void LyaNeutralHydrogenGasMix::performScattering(double lambda, const MaterialState* state, PhotonPacket* pp) const
{
    setScatteringInfoIfNeeded(const_cast<PhotonPacket*>(pp), state, lambda);
    auto scatinfo = const_cast<PhotonPacket*>(pp)->getScatteringInfo();

    // draw the outgoing direction from the dipole or the isotropic phase function
    // and, if required, update the polarization state of the photon packet
    Direction bfknew;
    if (scatinfo->dipole)
    {
        bfknew = _dpf.performScattering(pp->direction(), pp);
    }
    else
    {
        bfknew = random()->direction();
        if (includePolarization()) pp->setUnpolarized();
    }

    // Doppler-shift the photon packet wavelength into and out of the atom frame
    lambda = LyUtils::shiftWavelength(lambda, scatinfo->velocity, pp->direction(), bfknew);

    // execute the scattering event in the photon packet
    pp->scatter(bfknew, state->bulkVelocity(), lambda);
}

////////////////////////////////////////////////////////////////////

double LyaNeutralHydrogenGasMix::indicativeTemperature(const MaterialState* state, const Array& /*Jv*/) const
{
    return state->temperature();
}

////////////////////////////////////////////////////////////////////
