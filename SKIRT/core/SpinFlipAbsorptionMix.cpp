/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "SpinFlipAbsorptionMix.hpp"
#include "Configuration.hpp"
#include "Constants.hpp"
#include "DisjointWavelengthGrid.hpp"
#include "FatalError.hpp"
#include "MaterialState.hpp"

////////////////////////////////////////////////////////////////////

namespace
{
    // special wavelengths
    constexpr double lambdaSF = 21.10611405413e-2;  // 21 cm

    // wavelength range outside of which we consider absorption to be zero
    // (range of plus-min 9 dispersion elements using a velocity dispersion of 1000 km/s)
    constexpr Range absorptionRange(0.2047, 0.2174);

    // Einstein coefficient of the 21cm spin-flip transition
    constexpr double ASF = 2.8843e-15;
}

////////////////////////////////////////////////////////////////////

MaterialMix::MaterialType SpinFlipAbsorptionMix::materialType() const
{
    return MaterialType::Gas;
}

////////////////////////////////////////////////////////////////////

bool SpinFlipAbsorptionMix::hasExtraSpecificState() const
{
    return true;
}

////////////////////////////////////////////////////////////////////

vector<StateVariable> SpinFlipAbsorptionMix::specificStateVariableInfo() const
{
    return vector<StateVariable>{StateVariable::numberDensity(), StateVariable::temperature()};
}

////////////////////////////////////////////////////////////////////

void SpinFlipAbsorptionMix::initializeSpecificState(MaterialState* state, double /*metallicity*/, double temperature,
                                                    const Array& /*params*/) const
{
    // leave the properties at zero if the cell does not contain any material for this component
    if (state->numberDensity() > 0.)
    {
        // if no value was imported, use default value
        // make sure the temperature is at least the local universe CMB temperature
        state->setTemperature(max(Constants::Tcmb(), temperature >= 0. ? temperature : defaultTemperature()));
    }
}

////////////////////////////////////////////////////////////////////

double SpinFlipAbsorptionMix::mass() const
{
    return Constants::Mproton();
}

////////////////////////////////////////////////////////////////////

namespace
{
    // returns the absorption cross section per neutral hydrogen atom for the given wavelength and gas temperature
    double crossSection(double lambda, double T)
    {
        constexpr double front = 3. * M_SQRT2 * M_2_SQRTPI / M_PI / 128. * ASF * Constants::h() * Constants::c()
                                 * lambdaSF * lambdaSF / Constants::k();
        double Tspin = 6000. * (1 - exp(-0.0002 * T));
        double sigma = sqrt(Constants::k() / Constants::Mproton() * T);
        double u = Constants::c() * (lambda - lambdaSF) / lambda;
        double x = u / sigma;
        double expon = exp(-0.5 * x * x);
        return front / Tspin / sigma * expon;
    }
}

////////////////////////////////////////////////////////////////////

double SpinFlipAbsorptionMix::sectionAbs(double lambda) const
{
    return absorptionRange.contains(lambda) ? crossSection(lambda, defaultTemperature()) : 0.;
}

////////////////////////////////////////////////////////////////////

double SpinFlipAbsorptionMix::sectionSca(double /*lambda*/) const
{
    return 0.;
}

////////////////////////////////////////////////////////////////////

double SpinFlipAbsorptionMix::sectionExt(double lambda) const
{
    return sectionAbs(lambda);
}

////////////////////////////////////////////////////////////////////

double SpinFlipAbsorptionMix::opacityAbs(double lambda, const MaterialState* state, const PhotonPacket* /*pp*/) const
{
    if (absorptionRange.contains(lambda))
    {
        double number = state->numberDensity();
        if (number > 0.) return crossSection(lambda, state->temperature()) * number;
    }
    return 0.;
}

////////////////////////////////////////////////////////////////////

double SpinFlipAbsorptionMix::opacitySca(double /*lambda*/, const MaterialState* /*state*/,
                                         const PhotonPacket* /*pp*/) const
{
    return 0.;
}

////////////////////////////////////////////////////////////////////

double SpinFlipAbsorptionMix::opacityExt(double lambda, const MaterialState* state, const PhotonPacket* pp) const
{
    return opacityAbs(lambda, state, pp);
}

////////////////////////////////////////////////////////////////////

void SpinFlipAbsorptionMix::peeloffScattering(double& /*I*/, double& /*Q*/, double& /*U*/, double& /*V*/,
                                              double& /*lambda*/, Direction /*bfkobs*/, Direction /*bfky*/,
                                              const MaterialState* /*state*/, const PhotonPacket* /*pp*/) const
{}

////////////////////////////////////////////////////////////////////

void SpinFlipAbsorptionMix::performScattering(double /*lambda*/, const MaterialState* /*state*/,
                                              PhotonPacket* /*pp*/) const
{}

////////////////////////////////////////////////////////////////////

double SpinFlipAbsorptionMix::indicativeTemperature(const MaterialState* state, const Array& /*Jv*/) const
{
    return state->temperature();
}

////////////////////////////////////////////////////////////////////
