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
    // central wavelength
    constexpr double lambdaSF = Constants::lambdaSpinFlip();

    // wavelength range outside of which we consider absorption to be zero (approximately 20.47 - 21.74 cm)
    constexpr Range absorptionRange(lambdaSF * (1. - 0.03), lambdaSF * (1. + 0.03));
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
        constexpr double front = 3. * M_SQRT2 * M_2_SQRTPI / M_PI / 128. * Constants::EinsteinASpinFlip()
                                 * Constants::h() * Constants::c() * lambdaSF * lambdaSF / Constants::k();
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

double SpinFlipAbsorptionMix::indicativeTemperature(const MaterialState* state, const Array& /*Jv*/) const
{
    return state->temperature();
}

////////////////////////////////////////////////////////////////////
