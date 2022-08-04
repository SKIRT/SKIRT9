/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "SpinFlipHydrogenGasMix.hpp"
#include "Configuration.hpp"
#include "Constants.hpp"
#include "DisjointWavelengthGrid.hpp"
#include "FatalError.hpp"
#include "MaterialState.hpp"

////////////////////////////////////////////////////////////////////

namespace
{
    // special wavelengths
    constexpr double lambdaUV = 1000e-10;           // 1000 Angstrom
    constexpr double lambdaSF = 21.10611405413e-2;  // 21 cm

    // wavelength range outside of which we consider absorption to be zero (range of plus-min 0.21 mm)
    constexpr Range absorptionRange(lambdaSF * 0.999, lambdaSF * 1.001);

    // Einstein coefficient of the 21cm spin-flip transition
    constexpr double ASF = 2.8843e-15;

    // indices for custom state variables
    constexpr int NEUTRAL_FRACTION = 0;
    constexpr int ATOMIC_FRACTION = 1;
}

////////////////////////////////////////////////////////////////////

void SpinFlipHydrogenGasMix::setupSelfBefore()
{
    EmittingGasMix::setupSelfBefore();

    auto config = find<Configuration>();
    if (config->hasSecondaryEmission())
    {
        _indexUV = config->radiationFieldWLG()->bin(lambdaUV);
        if (_indexUV < 0) throw FATALERROR("Radiation field wavelength grid does not include 1000 Angstrom");
    }
}

////////////////////////////////////////////////////////////////////

bool SpinFlipHydrogenGasMix::hasExtraSpecificState() const
{
    return true;
}

////////////////////////////////////////////////////////////////////

MaterialMix::DynamicStateType SpinFlipHydrogenGasMix::hasDynamicMediumState() const
{
    return DynamicStateType::Secondary;
}

////////////////////////////////////////////////////////////////////

bool SpinFlipHydrogenGasMix::hasLineEmission() const
{
    return true;
}

////////////////////////////////////////////////////////////////////

vector<SnapshotParameter> SpinFlipHydrogenGasMix::parameterInfo() const
{
    return {SnapshotParameter::custom("neutral hydrogen fraction")};
}

////////////////////////////////////////////////////////////////////

vector<StateVariable> SpinFlipHydrogenGasMix::specificStateVariableInfo() const
{
    return vector<StateVariable>{StateVariable::numberDensity(), StateVariable::metallicity(),
                                 StateVariable::temperature(),
                                 StateVariable::custom(NEUTRAL_FRACTION, "neutral hydrogen fraction", ""),
                                 StateVariable::custom(ATOMIC_FRACTION, "atomic hydrogen fraction", "")};
}

////////////////////////////////////////////////////////////////////

void SpinFlipHydrogenGasMix::initializeSpecificState(MaterialState* state, double metallicity, double temperature,
                                                     const Array& params) const
{
    // leave the properties at zero if the cell does not contain any material for this component
    if (state->numberDensity() > 0.)
    {
        // if no value was imported, use default value
        // make sure the temperature is at least the local universe CMB temperature
        state->setMetallicity(metallicity >= 0. ? metallicity : defaultMetallicity());
        state->setTemperature(max(Constants::Tcmb(), temperature >= 0. ? temperature : defaultTemperature()));
        state->setCustom(NEUTRAL_FRACTION, params.size() ? params[0] : defaultNeutralFraction());
        state->setCustom(ATOMIC_FRACTION, 0.);
    }
}

////////////////////////////////////////////////////////////////////

UpdateStatus SpinFlipHydrogenGasMix::updateSpecificState(MaterialState* state, const Array& Jv) const
{
    if (_indexUV < 0) throw FATALERROR("State update should not be called if there is no radiation field");
    UpdateStatus status;

    // if the cell has no hydrogen or the neutral fraction is zero, then leave the atomic fraction at zero
    double nH = state->numberDensity();
    double fHIpH2 = state->custom(NEUTRAL_FRACTION);
    if (nH > 0. && fHIpH2 > 0.)
    {
        // get the radiation field
        // scaled to the reference Milky Way radiation field at 1000 Angstrom
        // converted from 1e6 photons/cm2/s/sr/eV to internal units W/m2/m/sr
        constexpr double h = Constants::h();
        constexpr double c = Constants::c();
        constexpr double Qel = Constants::Qelectron();
        constexpr double JMW = 1e6 * 1e4 * (h * c * h * c / (lambdaUV * lambdaUV * lambdaUV)) / Qel;
        double U = Jv[_indexUV] / JMW;

        // get the metallicity scaled to the solar reference value
        double D = state->metallicity() / 0.0127;

        // perform the partitioning scheme
        double Dstar = 1.5e-3 * log(1. + pow(3 * U, 1.7));
        double nstar = 25e6;
        double alpha = 2.5 * U / (1. + 0.25 * U * U);
        double s = 0.04 / (Dstar + D);
        double g = (1. + alpha * s + s * s) / (1. + s);
        double Lambda = log(1. + g * pow(D, 3. / 7.) * pow(U / 15., 4. / 7.));
        double x = pow(Lambda, 3. / 7.) * log(D / Lambda * nH / nstar);
        double fH2 = 1. / (1. + exp(-4. * x - 3. * x * x * x));

        // set the atomic fraction
        state->setCustom(ATOMIC_FRACTION, max(0., fHIpH2 - fH2));
        status.updateConverged();
    }
    return status;
}

////////////////////////////////////////////////////////////////////

double SpinFlipHydrogenGasMix::mass() const
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

double SpinFlipHydrogenGasMix::sectionAbs(double lambda) const
{
    return absorptionRange.contains(lambda) ? crossSection(lambda, defaultTemperature()) : 0.;
}

////////////////////////////////////////////////////////////////////

double SpinFlipHydrogenGasMix::sectionSca(double /*lambda*/) const
{
    return 0.;
}

////////////////////////////////////////////////////////////////////

double SpinFlipHydrogenGasMix::sectionExt(double lambda) const
{
    return sectionAbs(lambda);
}

////////////////////////////////////////////////////////////////////

double SpinFlipHydrogenGasMix::opacityAbs(double lambda, const MaterialState* state, const PhotonPacket* /*pp*/) const
{
    if (absorptionRange.contains(lambda))
    {
        double number = state->numberDensity() * state->custom(ATOMIC_FRACTION);
        if (number > 0.) return crossSection(lambda, state->temperature()) * number;
    }
    return 0.;
}

////////////////////////////////////////////////////////////////////

double SpinFlipHydrogenGasMix::opacitySca(double /*lambda*/, const MaterialState* /*state*/,
                                          const PhotonPacket* /*pp*/) const
{
    return 0.;
}

////////////////////////////////////////////////////////////////////

double SpinFlipHydrogenGasMix::opacityExt(double lambda, const MaterialState* state, const PhotonPacket* pp) const
{
    return opacityAbs(lambda, state, pp);
}

////////////////////////////////////////////////////////////////////

void SpinFlipHydrogenGasMix::peeloffScattering(double& /*I*/, double& /*Q*/, double& /*U*/, double& /*V*/,
                                               double& /*lambda*/, Direction /*bfkobs*/, Direction /*bfky*/,
                                               const MaterialState* /*state*/, const PhotonPacket* /*pp*/) const
{}

////////////////////////////////////////////////////////////////////

void SpinFlipHydrogenGasMix::performScattering(double /*lambda*/, const MaterialState* /*state*/,
                                               PhotonPacket* /*pp*/) const
{}

////////////////////////////////////////////////////////////////////

Array SpinFlipHydrogenGasMix::lineEmissionCenters() const
{
    Array centers(1);
    centers[0] = lambdaSF;
    return centers;
}

////////////////////////////////////////////////////////////////////

Array SpinFlipHydrogenGasMix::lineEmissionMasses() const
{
    Array masses(1);
    masses[0] = Constants::Mproton();
    return masses;
}

////////////////////////////////////////////////////////////////////

Array SpinFlipHydrogenGasMix::lineEmissionSpectrum(const MaterialState* state, const Array& /*Jv*/) const
{
    // calculate the 21 cm luminosity
    constexpr double front = 0.75 * ASF * Constants::h() * Constants::c() / lambdaSF;
    double L = front * state->custom(ATOMIC_FRACTION) * state->numberDensity() * state->volume();

    // encapsulate the result in an array
    Array luminosities(1);
    luminosities[0] = L;
    return luminosities;
}

////////////////////////////////////////////////////////////////////

double SpinFlipHydrogenGasMix::indicativeTemperature(const MaterialState* state, const Array& /*Jv*/) const
{
    return state->temperature();
}

////////////////////////////////////////////////////////////////////
