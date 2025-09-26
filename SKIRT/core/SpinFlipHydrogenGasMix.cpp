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
    constexpr double lambdaUV = 1000e-10;  // 1000 Angstrom
    constexpr double lambdaSF = Constants::lambdaSpinFlip();

    // wavelength range outside of which we consider absorption to be zero (approximately 20.47 - 21.74 cm)
    constexpr Range absorptionRange(lambdaSF * (1. - 0.03), lambdaSF * (1. + 0.03));

    // indices for custom state variables
    constexpr int NEUTRAL_SURFACE_DENSITY = 0;
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
    return {SnapshotParameter::custom("neutral hydrogen mass surface density", "masssurfacedensity", "Msun/pc2")};
}

////////////////////////////////////////////////////////////////////

vector<StateVariable> SpinFlipHydrogenGasMix::specificStateVariableInfo() const
{
    return vector<StateVariable>{
        StateVariable::numberDensity(), StateVariable::metallicity(), StateVariable::temperature(),
        StateVariable::custom(NEUTRAL_SURFACE_DENSITY, "neutral hydrogen mass surface density", "masssurfacedensity"),
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
        state->setCustom(NEUTRAL_SURFACE_DENSITY, params.size() ? params[0] : defaultNeutralSurfaceDensity());
        state->setCustom(ATOMIC_FRACTION, 0.);
    }
}

////////////////////////////////////////////////////////////////////

UpdateStatus SpinFlipHydrogenGasMix::updateSpecificState(MaterialState* state, const Array& Jv) const
{
    if (_indexUV < 0) throw FATALERROR("State update should not be called if there is no radiation field");
    UpdateStatus status;

    // if the cell has no neutral hydrogen, then leave the atomic fraction at zero
    if (state->numberDensity() > 0.)
    {
        // get the neutral hydrogen surface mass density
        double SigmaHIpH2 = state->custom(NEUTRAL_SURFACE_DENSITY);

        // get the dust-to-gas ratio (metallicity scaled to solar value)
        double DMW = state->metallicity() / 0.0127;

        // get the radiation field
        // scaled to the reference Milky Way radiation field at 1000 Angstrom
        // converted from 1e6 photons/cm2/s/sr/eV to internal units W/m2/m/sr
        constexpr double h = Constants::h();
        constexpr double c = Constants::c();
        constexpr double Qel = Constants::Qelectron();
        constexpr double JMW = 1e6 * 1e4 * (h * c * h * c / (lambdaUV * lambdaUV * lambdaUV)) / Qel;
        double UMW = Jv[_indexUV] / JMW;

        // get the length scale
        double Lcell = cbrt(state->volume());

        // perform the partitioning scheme
        constexpr double pc100 = 100. * Constants::pc();
        double S = Lcell / pc100;
        double S5 = pow(S, 5);
        double Dstar = 0.17 * (2. + S5) / (1. + S5);
        double g = sqrt(DMW * DMW + Dstar * Dstar);
        constexpr double front = 5e7 * Constants::Msun() / (1e6 * Constants::pc() * Constants::pc());
        double root = sqrt(0.001 + 0.1 * UMW);
        double Sigmac = front * root / (g * (1. + 1.69 * root));
        double alpha = 0.5 + 1. / (1. + sqrt(UMW * DMW * DMW / 600.));
        double Rmol = pow(SigmaHIpH2 / Sigmac, alpha);
        double fmol = Rmol / (Rmol + 1.);

        // set the atomic fraction
        state->setCustom(ATOMIC_FRACTION, max(0., 1. - fmol));
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
    constexpr double front = 0.75 * Constants::EinsteinASpinFlip() * Constants::h() * Constants::c() / lambdaSF;
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
