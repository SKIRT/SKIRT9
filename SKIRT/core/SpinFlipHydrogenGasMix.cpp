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
#include "PhotonPacket.hpp"

////////////////////////////////////////////////////////////////////

void SpinFlipHydrogenGasMix::setupSelfBefore()
{
    const double lambdaUV = 1000e-10;  // 1000 Angstrom

    auto config = find<Configuration>();
    if (config->hasSecondaryEmission())
    {
        _indexUV = config->radiationFieldWLG()->bin(lambdaUV);
        if (_indexUV < 0) throw FATALERROR("Radiation field wavelength grid does not include 1000 Angstrom");
    }
}

////////////////////////////////////////////////////////////////////

MaterialMix::MaterialType SpinFlipHydrogenGasMix::materialType() const
{
    return MaterialType::Gas;
}

////////////////////////////////////////////////////////////////////

bool SpinFlipHydrogenGasMix::hasExtraSpecificState() const
{
    return true;
}

////////////////////////////////////////////////////////////////////

bool SpinFlipHydrogenGasMix::hasSemiDynamicMediumState() const
{
    return true;
}

////////////////////////////////////////////////////////////////////

bool SpinFlipHydrogenGasMix::hasLineEmission() const
{
    return true;
}

////////////////////////////////////////////////////////////////////

vector<SnapshotParameter> SpinFlipHydrogenGasMix::parameterInfo() const
{
    return vector<SnapshotParameter>{SnapshotParameter("dust-to-gas ratio")};
}

////////////////////////////////////////////////////////////////////

vector<StateVariable> SpinFlipHydrogenGasMix::specificStateVariableInfo() const
{
    return vector<StateVariable>{StateVariable::numberDensity(), StateVariable::metallicity(),
                                 StateVariable::temperature(), StateVariable::custom(0, "dust-to-gas ratio", ""),
                                 StateVariable::custom(1, "UV field strength", "wavelengthmeanintensity")};
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
        state->setCustom(0, params.size() ? params[0] : defaultDustToGasRatio());
        state->setCustom(1, 0.);
    }
}

////////////////////////////////////////////////////////////////////

bool SpinFlipHydrogenGasMix::updateSpecificState(MaterialState* state, const Array& Jv) const
{
    if (_indexUV < 0) throw FATALERROR("State update should not be called if there is no radiation field");
    state->setCustom(1, Jv[_indexUV]);
    return true;
}

////////////////////////////////////////////////////////////////////

double SpinFlipHydrogenGasMix::mass() const
{
    return Constants::Mproton();
}

////////////////////////////////////////////////////////////////////

double SpinFlipHydrogenGasMix::sectionAbs(double /*lambda*/) const
{
    return 0.;  // TO DO
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

double SpinFlipHydrogenGasMix::opacityAbs(double /*lambda*/, const MaterialState* /*state*/,
                                          const PhotonPacket* /*pp*/) const
{
    return 0.;  // TO DO
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
                                               double& /*lambda*/, double /*w*/, Direction /*bfkobs*/,
                                               Direction /*bfky*/, const MaterialState* /*state*/,
                                               const PhotonPacket* /*pp*/) const
{}

////////////////////////////////////////////////////////////////////

void SpinFlipHydrogenGasMix::performScattering(double /*lambda*/, const MaterialState* /*state*/,
                                               PhotonPacket* /*pp*/) const
{}

////////////////////////////////////////////////////////////////////

Array SpinFlipHydrogenGasMix::lineEmissionCenters() const
{
    Array centers(1);
    centers[0] = 21.10611405413e-2;
    return centers;
}

////////////////////////////////////////////////////////////////////

double SpinFlipHydrogenGasMix::indicativeTemperature(const MaterialState* state, const Array& /*Jv*/) const
{
    return state->temperature();
}

////////////////////////////////////////////////////////////////////
