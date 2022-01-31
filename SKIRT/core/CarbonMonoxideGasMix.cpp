/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "CarbonMonoxideGasMix.hpp"
#include "Configuration.hpp"
#include "Constants.hpp"
#include "DisjointWavelengthGrid.hpp"
#include "FatalError.hpp"
#include "MaterialState.hpp"
#include "StringUtils.hpp"

////////////////////////////////////////////////////////////////////

namespace
{
    // number of level populations
    constexpr int numLevelPops = 40;

    // number of supported lines
    constexpr int numLines = 40;

    // ...
}

////////////////////////////////////////////////////////////////////

void CarbonMonoxideGasMix::setupSelfBefore()
{
    EmittingGasMix::setupSelfBefore();

    auto config = find<Configuration>();
    if (config->hasSecondaryEmission())
    {
        // ...
    }
}

////////////////////////////////////////////////////////////////////

bool CarbonMonoxideGasMix::hasExtraSpecificState() const
{
    return true;
}

////////////////////////////////////////////////////////////////////

bool CarbonMonoxideGasMix::hasSemiDynamicMediumState() const
{
    return true;
}

////////////////////////////////////////////////////////////////////

bool CarbonMonoxideGasMix::hasLineEmission() const
{
    return true;
}

////////////////////////////////////////////////////////////////////

vector<SnapshotParameter> CarbonMonoxideGasMix::parameterInfo() const
{
    return vector<SnapshotParameter>{SnapshotParameter("H2 number density", "volumenumberdensity", "1/cm3")};
}

////////////////////////////////////////////////////////////////////

vector<StateVariable> CarbonMonoxideGasMix::specificStateVariableInfo() const
{
    vector<StateVariable> result{StateVariable::numberDensity(), StateVariable::temperature(),
                                 StateVariable::custom(0, "H2 number density", "volumenumberdensity")};

    // add one custom variable for each level population
    // (the indices start at one because index zero is taken by the H2 number density)
    for (int p = 1; p <= numLevelPops; ++p)
        result.push_back(StateVariable::custom(p, "level population " + StringUtils::toString(p), ""));

    return result;
}

////////////////////////////////////////////////////////////////////

void CarbonMonoxideGasMix::initializeSpecificState(MaterialState* state, double /*metallicity*/, double temperature,
                                                   const Array& params) const
{
    // leave the properties untouched if the cell does not contain any material for this component
    if (state->numberDensity() > 0.)
    {
        // if no value was imported, use default value
        // make sure the temperature is at least the local universe CMB temperature
        state->setTemperature(max(Constants::Tcmb(), temperature >= 0. ? temperature : defaultTemperature()));
        state->setCustom(0, params.size() ? params[0] : state->numberDensity() * defaultMolecularHydrogenRatio());

        // clear the level populations
        for (int p = 1; p <= numLevelPops; ++p) state->setCustom(p, 0.);
    }
}

////////////////////////////////////////////////////////////////////

bool CarbonMonoxideGasMix::updateSpecificState(MaterialState* state, const Array& Jv) const
{
    // leave the properties untouched if the cell does not contain any material for this component
    if (state->numberDensity() > 0.)
    {
        // ...
        (void)Jv;

        // set the level populations (here some arbitrary example value
        for (int p = 1; p <= numLevelPops; ++p) state->setCustom(p, 1.);

        return true;
    }
    return false;
}

////////////////////////////////////////////////////////////////////

double CarbonMonoxideGasMix::mass() const
{
    return Constants::Mproton();  // ... should be CO molecule mass
}

////////////////////////////////////////////////////////////////////

double CarbonMonoxideGasMix::sectionAbs(double lambda) const
{
    // ...
    (void)lambda;

    return 1.;
}

////////////////////////////////////////////////////////////////////

double CarbonMonoxideGasMix::sectionSca(double /*lambda*/) const
{
    return 0.;
}

////////////////////////////////////////////////////////////////////

double CarbonMonoxideGasMix::sectionExt(double lambda) const
{
    return sectionAbs(lambda);
}

////////////////////////////////////////////////////////////////////

double CarbonMonoxideGasMix::opacityAbs(double lambda, const MaterialState* state, const PhotonPacket* /*pp*/) const
{
    // ...
    (void)lambda;
    (void)state;

    return 1.;
}

////////////////////////////////////////////////////////////////////

double CarbonMonoxideGasMix::opacitySca(double /*lambda*/, const MaterialState* /*state*/,
                                        const PhotonPacket* /*pp*/) const
{
    return 0.;
}

////////////////////////////////////////////////////////////////////

double CarbonMonoxideGasMix::opacityExt(double lambda, const MaterialState* state, const PhotonPacket* pp) const
{
    return opacityAbs(lambda, state, pp);
}

////////////////////////////////////////////////////////////////////

void CarbonMonoxideGasMix::peeloffScattering(double& /*I*/, double& /*Q*/, double& /*U*/, double& /*V*/,
                                             double& /*lambda*/, double /*w*/, Direction /*bfkobs*/, Direction /*bfky*/,
                                             const MaterialState* /*state*/, const PhotonPacket* /*pp*/) const
{}

////////////////////////////////////////////////////////////////////

void CarbonMonoxideGasMix::performScattering(double /*lambda*/, const MaterialState* /*state*/,
                                             PhotonPacket* /*pp*/) const
{}

////////////////////////////////////////////////////////////////////

Array CarbonMonoxideGasMix::lineEmissionCenters() const
{
    Array centers(numLines);
    // ...
    return centers;
}

////////////////////////////////////////////////////////////////////

Array CarbonMonoxideGasMix::lineEmissionMasses() const
{
    Array masses(numLines);
    // ...
    return masses;
}

////////////////////////////////////////////////////////////////////

Array CarbonMonoxideGasMix::lineEmissionSpectrum(const MaterialState* state, const Array& /*Jv*/) const
{
    Array luminosities(numLines);
    // ...
    (void)state;

    return luminosities;
}

////////////////////////////////////////////////////////////////////

double CarbonMonoxideGasMix::indicativeTemperature(const MaterialState* state, const Array& /*Jv*/) const
{
    return state->temperature();
}

////////////////////////////////////////////////////////////////////
