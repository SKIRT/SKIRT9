/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "SpinFlipHydrogenGasMix.hpp"
#include "Constants.hpp"
#include "MaterialState.hpp"
#include "PhotonPacket.hpp"

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

vector<SnapshotParameter> SpinFlipHydrogenGasMix::parameterInfo() const
{
    return vector<SnapshotParameter>{SnapshotParameter("dust-to-gas ratio")};
}

////////////////////////////////////////////////////////////////////

vector<StateVariable> SpinFlipHydrogenGasMix::specificStateVariableInfo() const
{
    return vector<StateVariable>{StateVariable::numberDensity(), StateVariable::metallicity(),
                                 StateVariable::temperature(), StateVariable::custom(0, "dust-to-gas ratio", "")};
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
        state->setMetallicity(metallicity > 0 ? metallicity : defaultMetallicity());
        state->setTemperature(max(Constants::Tcmb(), temperature > 0 ? temperature : defaultTemperature()));
        state->setCustom(0, params.size() && params[0] > 0 ? params[0] : defaultDustToGasRatio());
    }
}

////////////////////////////////////////////////////////////////////

double SpinFlipHydrogenGasMix::mass() const
{
    return Constants::Mproton();
}

////////////////////////////////////////////////////////////////////

double SpinFlipHydrogenGasMix::sectionAbs(double /*lambda*/) const
{
    return 0.;
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
                                               double& /*lambda*/, double /*w*/, Direction /*bfkobs*/,
                                               Direction /*bfky*/, const MaterialState* /*state*/,
                                               const PhotonPacket* /*pp*/) const
{}

////////////////////////////////////////////////////////////////////

void SpinFlipHydrogenGasMix::performScattering(double /*lambda*/, const MaterialState* /*state*/,
                                               PhotonPacket* /*pp*/) const
{}

////////////////////////////////////////////////////////////////////

double SpinFlipHydrogenGasMix::indicativeTemperature(const MaterialState* state, const Array& /*Jv*/) const
{
    return state->temperature();
}

////////////////////////////////////////////////////////////////////
