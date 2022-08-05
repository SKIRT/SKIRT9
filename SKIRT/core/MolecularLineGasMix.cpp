/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "MolecularLineGasMix.hpp"
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

void MolecularLineGasMix::setupSelfBefore()
{
    EmittingGasMix::setupSelfBefore();

    auto config = find<Configuration>();
    if (config->hasSecondaryEmission())
    {
        // ...
    }
}

////////////////////////////////////////////////////////////////////

bool MolecularLineGasMix::hasNegativeExtinction() const
{
    return true;
}

////////////////////////////////////////////////////////////////////

bool MolecularLineGasMix::hasExtraSpecificState() const
{
    return true;
}

////////////////////////////////////////////////////////////////////

MaterialMix::DynamicStateType MolecularLineGasMix::hasDynamicMediumState() const
{
    return DynamicStateType::PrimaryIfMergedIterations;
}

////////////////////////////////////////////////////////////////////

bool MolecularLineGasMix::hasLineEmission() const
{
    return true;
}

////////////////////////////////////////////////////////////////////

vector<SnapshotParameter> MolecularLineGasMix::parameterInfo() const
{
    return {SnapshotParameter::custom("H2 number density", "numbervolumedensity", "1/cm3")};
}

////////////////////////////////////////////////////////////////////

vector<StateVariable> MolecularLineGasMix::specificStateVariableInfo() const
{
    vector<StateVariable> result{StateVariable::numberDensity(), StateVariable::temperature(),
                                 StateVariable::custom(0, "H2 number density", "numbervolumedensity")};

    // add one custom variable for each level population
    // (the indices start at one because index zero is taken by the H2 number density)
    for (int p = 1; p <= numLevelPops; ++p)
        result.push_back(StateVariable::custom(p, "level population " + StringUtils::toString(p), ""));

    return result;
}

////////////////////////////////////////////////////////////////////

void MolecularLineGasMix::initializeSpecificState(MaterialState* state, double /*metallicity*/, double temperature,
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

UpdateStatus MolecularLineGasMix::updateSpecificState(MaterialState* state, const Array& Jv) const
{
    UpdateStatus status;

    // leave the properties untouched if the cell does not contain any material for this component
    if (state->numberDensity() > 0.)
    {
        // ...
        (void)Jv;

        // set the level populations (here some arbitrary example value
        for (int p = 1; p <= numLevelPops; ++p) state->setCustom(p, 1.);

        status.updateConverged();
    }
    return status;
}

////////////////////////////////////////////////////////////////////

bool MolecularLineGasMix::isSpecificStateConverged(int numCells, int /*numUpdated*/, int numNotConverged) const
{
    return static_cast<double>(numNotConverged) / static_cast<double>(numCells) <= maxFractionNotConverged();
}

////////////////////////////////////////////////////////////////////

double MolecularLineGasMix::mass() const
{
    return Constants::Mproton();  // ... should be the mass of the species
}

////////////////////////////////////////////////////////////////////

double MolecularLineGasMix::sectionAbs(double /*lambda*/) const
{
    return 0.;
}

////////////////////////////////////////////////////////////////////

double MolecularLineGasMix::sectionSca(double /*lambda*/) const
{
    return 0.;
}

////////////////////////////////////////////////////////////////////

double MolecularLineGasMix::sectionExt(double /*lambda*/) const
{
    return 0.;
}

////////////////////////////////////////////////////////////////////

double MolecularLineGasMix::opacityAbs(double lambda, const MaterialState* state, const PhotonPacket* /*pp*/) const
{
    // ...
    (void)lambda;
    (void)state;

    return 1.;
}

////////////////////////////////////////////////////////////////////

double MolecularLineGasMix::opacitySca(double /*lambda*/, const MaterialState* /*state*/,
                                       const PhotonPacket* /*pp*/) const
{
    return 0.;
}

////////////////////////////////////////////////////////////////////

double MolecularLineGasMix::opacityExt(double lambda, const MaterialState* state, const PhotonPacket* pp) const
{
    return opacityAbs(lambda, state, pp);
}

////////////////////////////////////////////////////////////////////

void MolecularLineGasMix::peeloffScattering(double& /*I*/, double& /*Q*/, double& /*U*/, double& /*V*/,
                                            double& /*lambda*/, Direction /*bfkobs*/, Direction /*bfky*/,
                                            const MaterialState* /*state*/, const PhotonPacket* /*pp*/) const
{}

////////////////////////////////////////////////////////////////////

void MolecularLineGasMix::performScattering(double /*lambda*/, const MaterialState* /*state*/,
                                            PhotonPacket* /*pp*/) const
{}

////////////////////////////////////////////////////////////////////

Array MolecularLineGasMix::lineEmissionCenters() const
{
    Array centers(numLines);
    for (int i = 0; i != numLines; ++i) centers[i] = i + 1;
    // ...
    return centers;
}

////////////////////////////////////////////////////////////////////

Array MolecularLineGasMix::lineEmissionMasses() const
{
    Array masses(numLines);
    // ...
    return masses;
}

////////////////////////////////////////////////////////////////////

Array MolecularLineGasMix::lineEmissionSpectrum(const MaterialState* state, const Array& /*Jv*/) const
{
    Array luminosities(numLines);
    // ...
    (void)state;

    return luminosities;
}

////////////////////////////////////////////////////////////////////

double MolecularLineGasMix::indicativeTemperature(const MaterialState* state, const Array& /*Jv*/) const
{
    return state->temperature();
}

////////////////////////////////////////////////////////////////////
