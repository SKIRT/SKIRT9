/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "NonLTELineGasMix.hpp"
#include "Configuration.hpp"
#include "Constants.hpp"
#include "DisjointWavelengthGrid.hpp"
#include "FatalError.hpp"
#include "Log.hpp"
#include "MaterialState.hpp"
#include "NR.hpp"
#include "StringUtils.hpp"
#include "TextInFile.hpp"
#include "Units.hpp"

////////////////////////////////////////////////////////////////////

void NonLTELineGasMix::setupSelfBefore()
{
    EmittingGasMix::setupSelfBefore();

    // get the name of the configured species and the name(s) of its collision partners
    string name;
    vector<string> colNames{"H2"};  // just molecular hydrogen by default
    switch (species())
    {
        // molecular species
        case Species::Test: name = "TT"; break;
        case Species::Hydroxyl:
            name = "OH";
            colNames = {"H2", "H", "He"};
            break;
        case Species::HydroxylHFS: name = "OHhfs"; break;
        case Species::Formyl: name = "HCO+"; break;
        case Species::HydrogenCyanide:
            name = "HCN";
            colNames = {"H2", "e-"};
            break;
        case Species::CarbonMonoxide: name = "CO"; break;
        case Species::MolecularHydrogen:
            name = "H2";
            colNames = {"H2", "H", "H+", "He"};
            break;

        // atomic species
        case Species::AtomicCarbon:
            name = "C_I";
            colNames = {"H2", "H", "H+", "e-", "He"};
            break;
        case Species::IonizedCarbon:
            name = "C_II";
            colNames = {"H2", "H", "e-"};
            break;
        case Species::DoublyIonizedCarbon:
            name = "C_III";
            colNames = {"e-"};
            break;
        case Species::TriplyIonizedCarbon:
            name = "C_IV";
            colNames = {"e-"};
            break;
        case Species::FourTimesIonizedCarbon:
            name = "C_V";
            colNames = {"e-"};
            break;
        case Species::FiveTimesIonizedCarbon:
            name = "C_VI";
            colNames = {"e-"};
            break;

        case Species::AtomicNitrogen:
            name = "N_I";
            colNames = {"e-"};
            break;
        case Species::IonizedNitrogen:
            name = "N_II";
            colNames = {"H", "e-"};
            break;
        case Species::DoublyIonizedNitrogen:
            name = "N_III";
            colNames = {"e-"};
            break;
        case Species::TriplyIonizedNitrogen:
            name = "N_IV";
            colNames = {"e-"};
            break;
        case Species::FourTimesIonizedNitrogen:
            name = "N_V";
            colNames = {"e-"};
            break;
        case Species::FiveTimesIonizedNitrogen:
            name = "N_VI";
            colNames = {"e-"};
            break;
        case Species::SixTimesIonizedNitrogen:
            name = "N_VII";
            colNames = {"e-"};
            break;

        case Species::AtomicOxygen:
            name = "O_I";
            colNames = {"H2", "H", "H+", "e-", "He"};
            break;
        case Species::IonizedOxygen:
            name = "O_II";
            colNames = {"e-"};
            break;
        case Species::DoublyIonizedOxygen:
            name = "O_III";
            colNames = {"H", "e-"};
            break;
        case Species::TriplyIonizedOxygen:
            name = "O_IV";
            colNames = {"e-"};
            break;
        case Species::FourTimesIonizedOxygen:
            name = "O_V";
            colNames = {"e-"};
            break;
        case Species::FiveTimesIonizedOxygen:
            name = "O_VI";
            colNames = {"e-"};
            break;
        case Species::SixTimesIonizedOxygen:
            name = "O_VII";
            colNames = {"e-"};
            break;
        case Species::SevenTimesIonizedOxygen:
            name = "O_VIII";
            colNames = {"e-"};
            break;

        case Species::AtomicNeon:
            name = "Ne_I";
            colNames = {"e-"};
            break;
        case Species::IonizedNeon:
            name = "Ne_II";
            colNames = {"e-"};
            break;
        case Species::DoublyIonizedNeon:
            name = "Ne_III";
            colNames = {"e-"};
            break;
        case Species::TriplyIonizedNeon:
            name = "Ne_IV";
            colNames = {"e-"};
            break;
        case Species::FourTimesIonizedNeon:
            name = "Ne_V";
            colNames = {"e-"};
            break;
        case Species::FiveTimesIonizedNeon:
            name = "Ne_VI";
            colNames = {"e-"};
            break;
        case Species::SixTimesIonizedNeon:
            name = "Ne_VII";
            colNames = {"e-"};
            break;
        case Species::SevenTimesIonizedNeon:
            name = "Ne_VIII";
            colNames = {"e-"};
            break;
        case Species::EightTimesIonizedNeon:
            name = "Ne_IX";
            colNames = {"e-"};
            break;
        case Species::NineTimesIonizedNeon:
            name = "Ne_X";
            colNames = {"e-"};
            break;

        case Species::IonizedSodium:
            name = "Na_II";
            colNames = {"e-"};
            break;
        case Species::DoublyIonizedSodium:
            name = "Na_III";
            colNames = {"e-"};
            break;
        case Species::TriplyIonizedSodium:
            name = "Na_IV";
            colNames = {"e-"};
            break;
        case Species::FourTimesIonizedSodium:
            name = "Na_V";
            colNames = {"e-"};
            break;
        case Species::FiveTimesIonizedSodium:
            name = "Na_VI";
            colNames = {"e-"};
            break;
        case Species::SixTimesIonizedSodium:
            name = "Na_VII";
            colNames = {"e-"};
            break;
        case Species::SevenTimesIonizedSodium:
            name = "Na_VIII";
            colNames = {"e-"};
            break;
        case Species::EightTimesIonizedSodium:
            name = "Na_IX";
            colNames = {"e-"};
            break;
        case Species::NineTimesIonizedSodium:
            name = "Na_X";
            colNames = {"e-"};
            break;
        case Species::TenTimesIonizedSodium:
            name = "Na_XI";
            colNames = {"e-"};
            break;

        case Species::IonizedMagnesium:
            name = "Mg_II";
            colNames = {"e-"};
            break;
        case Species::DoublyIonizedMagnesium:
            name = "Mg_III";
            colNames = {"e-"};
            break;
        case Species::TriplyIonizedMagnesium:
            name = "Mg_IV";
            colNames = {"e-"};
            break;
        case Species::FourTimesIonizedMagnesium:
            name = "Mg_V";
            colNames = {"e-"};
            break;
        case Species::FiveTimesIonizedMagnesium:
            name = "Mg_VI";
            colNames = {"e-"};
            break;
        case Species::SixTimesIonizedMagnesium:
            name = "Mg_VII";
            colNames = {"e-"};
            break;
        case Species::SevenTimesIonizedMagnesium:
            name = "Mg_VIII";
            colNames = {"e-"};
            break;
        case Species::EightTimesIonizedMagnesium:
            name = "Mg_IX";
            colNames = {"e-"};
            break;
        case Species::NineTimesIonizedMagnesium:
            name = "Mg_X";
            colNames = {"e-"};
            break;
        case Species::TenTimesIonizedMagnesium:
            name = "Mg_XI";
            colNames = {"e-"};
            break;

        case Species::AtomicSilicon:
            name = "Si_I";
            colNames = {"H", "e-"};
            break;
        case Species::IonizedSilicon:
            name = "Si_II";
            colNames = {"e-"};
            break;
        case Species::DoublyIonizedSilicon:
            name = "Si_III";
            colNames = {"e-"};
            break;
        case Species::TriplyIonizedSilicon:
            name = "Si_IV";
            colNames = {"e-"};
            break;
        case Species::FourTimesIonizedSilicon:
            name = "Si_V";
            colNames = {"e-"};
            break;
        case Species::FiveTimesIonizedSilicon:
            name = "Si_VI";
            colNames = {"e-"};
            break;
        case Species::SixTimesIonizedSilicon:
            name = "Si_VII";
            colNames = {"e-"};
            break;
        case Species::SevenTimesIonizedSilicon:
            name = "Si_VIII";
            colNames = {"e-"};
            break;
        case Species::EightTimesIonizedSilicon:
            name = "Si_IX";
            colNames = {"e-"};
            break;
        case Species::NineTimesIonizedSilicon:
            name = "Si_X";
            colNames = {"e-"};
            break;
        case Species::TenTimesIonizedSilicon:
            name = "Si_XI";
            colNames = {"e-"};
            break;

        case Species::AtomicSulfur:
            name = "S_I";
            colNames = {"H", "e-"};
            break;
        case Species::IonizedSulfur:
            name = "S_II";
            colNames = {"e-"};
            break;
        case Species::DoublyIonizedSulfur:
            name = "S_III";
            colNames = {"e-"};
            break;
        case Species::TriplyIonizedSulfur:
            name = "S_IV";
            colNames = {"e-"};
            break;
        case Species::FourTimesIonizedSulfur:
            name = "S_V";
            colNames = {"e-"};
            break;
        case Species::FiveTimesIonizedSulfur:
            name = "S_VI";
            colNames = {"e-"};
            break;
        case Species::SixTimesIonizedSulfur:
            name = "S_VII";
            colNames = {"e-"};
            break;
        case Species::SevenTimesIonizedSulfur:
            name = "S_VIII";
            colNames = {"e-"};
            break;
        case Species::EightTimesIonizedSulfur:
            name = "S_IX";
            colNames = {"e-"};
            break;
        case Species::NineTimesIonizedSulfur:
            name = "S_X";
            colNames = {"e-"};
            break;
        case Species::TenTimesIonizedSulfur:
            name = "S_XI";
            colNames = {"e-"};
            break;

        case Species::IonizedIron:
            name = "Fe_II";
            colNames = {"e-"};
            break;
        case Species::DoublyIonizedIron:
            name = "Fe_III";
            colNames = {"e-"};
            break;
        case Species::TriplyIonizedIron:
            name = "Fe_IV";
            colNames = {"e-"};
            break;
        case Species::FourTimesIonizedIron:
            name = "Fe_V";
            colNames = {"e-"};
            break;
        case Species::FiveTimesIonizedIron:
            name = "Fe_VI";
            colNames = {"e-"};
            break;
        case Species::SixTimesIonizedIron:
            name = "Fe_VII";
            colNames = {"e-"};
            break;
        case Species::SevenTimesIonizedIron:
            name = "Fe_VIII";
            colNames = {"e-"};
            break;
        case Species::EightTimesIonizedIron:
            name = "Fe_IX";
            colNames = {"e-"};
            break;
        case Species::NineTimesIonizedIron:
            name = "Fe_X";
            colNames = {"e-"};
            break;
        case Species::TenTimesIonizedIron:
            name = "Fe_XI";
            colNames = {"e-"};
            break;

        case Species::DoublyIonizedArgon:
            name = "Ar_III";
            colNames = {"e-"};
            break;
        case Species::TriplyIonizedArgon:
            name = "Ar_IV";
            colNames = {"e-"};
            break;
        case Species::FourTimesIonizedArgon:
            name = "Ar_V";
            colNames = {"e-"};
            break;
        case Species::FiveTimesIonizedArgon:
            name = "Ar_VI";
            colNames = {"e-"};
            break;
        case Species::SixTimesIonizedArgon:
            name = "Ar_VII";
            colNames = {"e-"};
            break;
        case Species::SevenTimesIonizedArgon:
            name = "Ar_VIII";
            colNames = {"e-"};
            break;
        case Species::EightTimesIonizedArgon:
            name = "Ar_IX";
            colNames = {"e-"};
            break;
        case Species::NineTimesIonizedArgon:
            name = "Ar_X";
            colNames = {"e-"};
            break;
        case Species::TenTimesIonizedArgon:
            name = "Ar_XI";
            colNames = {"e-"};
            break;
    }
    _name = name;

    // load the atomic model (mass, energy levels, radiative and collisional transitions) using the
    // solver's shared loader; this issues no per-file log messages of its own (see loadAtomicModel()),
    // so a single summary message is issued below instead
    _gasLineEmission.initialize(this);
    _gasLineEmission.loadAtomicModel(name, colNames, numEnergyLevels(), _model);
    auto log = find<Log>();
    log->info("Loaded atomic model for " + name + " from " + std::to_string(3 + 2 * colNames.size())
              + " resource files");

    // log summary info on the radiative lines
    auto units = find<Units>();
    int numLines = _model.numLines();
    log->info("Radiative lines for " + name + ":");
    if (numLines == 0) throw FATALERROR("There are no radiative transitions; increase the number of energy levels");
    for (int k = 0; k != numLines; ++k)
    {
        if (_model.branchRatio[k] > lowestBranchingRatio() && _model.indexUpRad[k] < maxUpperLevelForRadiation())
        {
            log->info("  (" + StringUtils::toString(_model.indexUpRad[k]) + "-"
                      + StringUtils::toString(_model.indexLowRad[k]) + ") "
                      + StringUtils::toString(units->owavelength(_model.center[k])) + " " + units->uwavelength()
                      + ", branch ratio=" + StringUtils::toString(_model.branchRatio[k]) + ">"
                      + StringUtils::toString(lowestBranchingRatio()) + ": radiative transition included");
        }
        else
        {
            log->info("  (" + StringUtils::toString(_model.indexUpRad[k]) + "-"
                      + StringUtils::toString(_model.indexLowRad[k]) + ") "
                      + StringUtils::toString(units->owavelength(_model.center[k])) + " " + units->uwavelength()
                      + ", branch ratio=" + StringUtils::toString(_model.branchRatio[k]) + "<"
                      + StringUtils::toString(lowestBranchingRatio()) + ": radiative transition excluded");
        }
    }

    // log summary info on the collisional partner(s)
    log->info("Collisional partner(s) for " + name + ": " + StringUtils::join(colNames, ", "));

    // verify that the radiation field wavelength grid, if present, has a bin covering the line centers
    // and cache the characteristic wavelengths and bin widths
    auto rfwlg = find<Configuration>()->radiationFieldWLG();
    if (rfwlg)
    {
        rfwlg->setup();
        for (int k = 0; k != numLines; ++k)
        {
            if (_model.branchRatio[k] > lowestBranchingRatio() && _model.indexUpRad[k] < maxUpperLevelForRadiation())
            {
                if (rfwlg->bin(_model.center[k]) < 0)
                    throw FATALERROR("Radiation field wavelength grid does not cover the central line for transition ("
                                     + StringUtils::toString(_model.indexUpRad[k]) + "-"
                                     + StringUtils::toString(_model.indexLowRad[k]) + ")");
            }
        }
        _numWavelengths = rfwlg->numBins();
        _lambdav = rfwlg->lambdav();
        _dlambdav = rfwlg->dlambdav();
    }

    // load the initial relative level populations if the user provided a filename
    if (!initialLevelPopsFilename().empty())
    {
        TextInFile infile(this, initialLevelPopsFilename(), "initial level populations");
        infile.addColumn("cell index");
        for (int p = 0; p != _model.numLevels(); ++p)
            infile.addColumn("population of level " + std::to_string(p), "numbervolumedensity", "1/cm3");
        _initLevelPops = infile.readAllRows();
    }
}

////////////////////////////////////////////////////////////////////

bool NonLTELineGasMix::hasNegativeExtinction() const
{
    return true;
}

////////////////////////////////////////////////////////////////////

bool NonLTELineGasMix::hasExtraSpecificState() const
{
    return true;
}

////////////////////////////////////////////////////////////////////

MaterialMix::DynamicStateType NonLTELineGasMix::hasDynamicMediumState() const
{
    return DynamicStateType::PrimaryIfMergedIterations;
}

////////////////////////////////////////////////////////////////////

bool NonLTELineGasMix::hasLineEmission() const
{
    return true;
}

////////////////////////////////////////////////////////////////////

vector<SnapshotParameter> NonLTELineGasMix::parameterInfo() const
{
    vector<SnapshotParameter> result;

    // add the number density of each collisional partner
    for (const auto& partner : _model.colPartner)
        result.push_back(SnapshotParameter::custom(partner.name + " number density", "numbervolumedensity", "1/cm3"));

    // add the turbulence velocity
    result.push_back(SnapshotParameter::custom("turbulence velocity", "velocity", "km/s"));

    return result;
}

////////////////////////////////////////////////////////////////////

vector<StateVariable> NonLTELineGasMix::specificStateVariableInfo() const
{
    // add standard variables for the number density of the species under consideration
    // and for the effective gas temperature (including kinetic temperature and unresolved turbulence)
    vector<StateVariable> result{StateVariable::numberDensity(), StateVariable::temperature()};

    // next available custom variable index
    int index = 0;

    // add custom variable for the kinetic gas temperature (i.e. excluding turbulence)
    const_cast<NonLTELineGasMix*>(this)->_indexKineticTemperature = index;
    result.push_back(StateVariable::custom(index++, "kinetic gas temperature", "temperature"));

    // add custom variable for the number density of each collisional partner
    const_cast<NonLTELineGasMix*>(this)->_indexFirstColPartnerDensity = index;
    for (const auto& partner : _model.colPartner)
        result.push_back(StateVariable::custom(index++, partner.name + " number density", "numbervolumedensity"));

    // add custom variable for the population of each energy level
    const_cast<NonLTELineGasMix*>(this)->_indexFirstLevelPopulation = index;
    for (int p = 0; p != _model.numLevels(); ++p)
        result.push_back(
            StateVariable::custom(index++, "population of level " + std::to_string(p), "numbervolumedensity"));

    // if requested, add custom variable for the line-profile-averaged mean intensity at each transition line
    if (storeMeanIntensities())
    {
        const_cast<NonLTELineGasMix*>(this)->_indexFirstMeanIntensity = index;
        for (int k = 0; k != _model.numLines(); ++k)
            result.push_back(StateVariable::custom(index++, "mean intensity at line " + std::to_string(k),
                                                   "wavelengthmeanintensity"));
    }
    return result;
}

////////////////////////////////////////////////////////////////////

// Macro's for accessing custom variables in the material state
// This is an ugly hack but there does not seem to be an elegant in-language mechanism to accomplish this
#define setKineticTemperature(value) setCustom(_indexKineticTemperature, (value))
#define kineticTemperature() custom(_indexKineticTemperature)
#define setColPartnerDensity(index, value) setCustom(_indexFirstColPartnerDensity + (index), (value))
#define colPartnerDensity(index) custom(_indexFirstColPartnerDensity + (index))
#define setLevelPopulation(index, value) setCustom(_indexFirstLevelPopulation + (index), (value))
#define levelPopulation(index) custom(_indexFirstLevelPopulation + (index))
#define setMeanIntensity(index, value) setCustom(_indexFirstMeanIntensity + (index), (value))

////////////////////////////////////////////////////////////////////

void NonLTELineGasMix::initializeSpecificState(MaterialState* state, double /*metallicity*/, double temperature,
                                               const Array& params) const
{
    // if the cell does not contain any material for this component, leave all properties at zero values
    if (state->numberDensity() > 0.)
    {
        int numColPartners = _model.numColPartners();
        int numLevels = _model.numLevels();

        // copy kinetic temperature from import or default
        double Tkin = temperature >= 0. ? temperature : defaultTemperature();
        state->setKineticTemperature(Tkin);

        // set effective temperature, including imported or default turbulence
        double vturb = params.size() ? params[numColPartners] : defaultTurbulenceVelocity();
        double Teff = Tkin + 0.5 * vturb * vturb * _model.mass / Constants::k();
        state->setTemperature(Teff);

        // copy collisional partner densities from import or default
        if (params.size())
        {
            for (int c = 0; c != numColPartners; ++c) state->setColPartnerDensity(c, params[c]);
        }
        else
        {
            const auto& ratios = defaultCollisionPartnerRatios();
            if (static_cast<int>(ratios.size()) < numColPartners)
                throw FATALERROR("The number of collision partners exceeds the number of default ratios");
            for (int c = 0; c != numColPartners; ++c)
                state->setColPartnerDensity(c, state->numberDensity() * ratios[c]);
        }

        // initialize level population using boltzmann distribution (i.e., start with LTE)
        auto initLTE = [this, Tkin, state, numLevels]() {
            Array levelPops(numLevels);
            for (int p = 0; p != numLevels; ++p)
                levelPops[p] = _model.weight[p] * exp(-_model.energy[p] / Constants::k() / Tkin);
            // normalize and store; a non-positive Tkin (schema-legal: defaultTemperature allows 0, and
            // so does an imported temperature) makes the sum zero or NaN, so fall back to putting all
            // population in the ground state -- the T->0 limit of a Boltzmann distribution -- instead
            double sum = levelPops.sum();
            if (sum > 0.)
            {
                levelPops *= state->numberDensity() / sum;
            }
            else
            {
                for (int p = 0; p != numLevels; ++p) levelPops[p] = 0.;
                levelPops[0] = state->numberDensity();
            }
            for (int p = 0; p != numLevels; ++p) state->setLevelPopulation(p, levelPops[p]);
        };

        if (initialLevelPopsCase() == InitialLevelPopsCase::LTE)
        {
            initLTE();
        }
        else if (initialLevelPopsCase() == InitialLevelPopsCase::CollisionallyExcited)
        {
            // initilize level population based on the the detailed balance condition
            // between collisional transitions and spontaneous emission (i.e., start with non-LTE);
            // call the solver directly (rather than through updateSpecificState()) so that this
            // initialization always happens, even if updateDynamicStatesFlag disables later updates
            Array Jv_zero(_numWavelengths);  // automatically initialized to zero values
            solveLevelPopulations(state, Jv_zero);
        }
        else if (initialLevelPopsCase() == InitialLevelPopsCase::Custom)
        {
            // if the user configured a file with initial level populations, use those data instead;
            // fall back to LTE for any cell whose index is not present in that file, or whose row
            // does not sum to a positive value
            size_t m = state->cellIndex();
            Array levelPops(numLevels);
            double sum = 0.;
            if (m < _initLevelPops.size())
            {
                for (int p = 0; p != numLevels; ++p) levelPops[p] = _initLevelPops[m][p + 1];
                sum = levelPops.sum();
            }
            if (sum > 0.)
            {
                // normalize and store
                levelPops *= state->numberDensity() / sum;
                for (int p = 0; p != numLevels; ++p) state->setLevelPopulation(p, levelPops[p]);
            }
            else
            {
                initLTE();
            }
        }
    }
}

////////////////////////////////////////////////////////////////////

namespace
{
    // return the dispersion of a line profile in wavelength space
    // given the line center, the effective gas temperature, and the species mass
    double sigmaForLine(double center, double temperature, double mass)
    {
        return center / Constants::c() * sqrt(Constants::k() * temperature / mass);
    }

    // return the value at ordinate x of a normalized Gaussian probability distribution
    // with given center mu and dispersion sigma
    double gaussian(double x, double mu, double sigma)
    {
        double u = (x - mu) / sigma;
        constexpr double front = 0.25 * M_SQRT2 * M_2_SQRTPI;
        return front / sigma * exp(-0.5 * u * u);
    }

    // hardcoded constant indicating the line profile range considered in the calculations
    // expressed as a multiple of the Gaussian sigma (in each direction from the center)
    constexpr double PROFILE_RANGE = 4.;

    // hardcoded constant indicating the fractional error allowed on the integration of
    // the Gaussian line profile over the simulation's radiation field wavelength grid
    constexpr double MAX_GAUSS_ERROR_WARN = 0.01;
    constexpr double MAX_GAUSS_ERROR_FAIL = 0.10;
}

////////////////////////////////////////////////////////////////////

UpdateStatus NonLTELineGasMix::updateSpecificState(MaterialState* state, const Array& Jv) const
{
    // initialize status indicator to "not updated"
    UpdateStatus status;

    // if the cell does not contain any material for this component, leave all properties untouched
    if (updateDynamicStatesFlag() == false)
    {
        status.updateConverged();
    }
    else if (state->numberDensity() > 0)
    {
        double change = solveLevelPopulations(state, Jv);

        // verify convergence
        if (change > maxChangeInLevelPopulations())
            status.updateNotConverged();
        else
            status.updateConverged();
    }
    return status;
}

////////////////////////////////////////////////////////////////////

double NonLTELineGasMix::solveLevelPopulations(MaterialState* state, const Array& Jv) const
{
    // gather the per-cell inputs for the shared level-population solver; the Einstein A
    // (spontaneous emission) term is always included by the shared solver, but the stimulated
    // Bul/Blu terms only act through meanJ, so a line excluded below simply keeps meanJ at its
    // default of zero rather than being removed from the model
    GasLineEmission::Environment env;
    env.Tkin = state->kineticTemperature();
    env.nTotal = state->numberDensity();
    int numColPartners = _model.numColPartners();
    int numLines = _model.numLines();
    int numLevels = _model.numLevels();
    env.nPartner.resize(numColPartners);
    for (int c = 0; c != numColPartners; ++c) env.nPartner[c] = state->colPartnerDensity(c);
    env.meanJ.assign(numLines, 0.);

    // calculate the mean intensity of the radiation field convolved over the normalized line profile g
    // for each radiative transition:
    //   J_convolved = \int J_lambda(lambda) g(lambda) d lambda  /  \int g(lambda) d lambda
    // all wavelength points within a given range around the line center are used, and the grid is
    // verified to be fine enough to reproduce the normalization 1 = \int g(lambda) d lambda
    for (int k = 0; k != numLines; ++k)
    {
        int up = _model.indexUpRad[k];
        int low = _model.indexLowRad[k];

        // ignore radiative transitions from high upper levels.
        if (up >= maxUpperLevelForRadiation()) continue;

        // ignore radiative transitions with a branching ratio below the threshold
        if (_model.branchRatio[k] < lowestBranchingRatio()) continue;

        auto log = find<Log>();
        double center = _model.center[k];
        double sigma = sigmaForLine(center, state->temperature(), _model.mass);
        double lambdamin = center - PROFILE_RANGE * sigma;
        double lambdamax = center + PROFILE_RANGE * sigma;
        int ellmin = std::lower_bound(begin(_lambdav), end(_lambdav), lambdamin) - begin(_lambdav);
        int ellmax = std::upper_bound(begin(_lambdav), end(_lambdav), lambdamax) - begin(_lambdav);
        double gsum = 0.;
        double Jsum = 0.;
        for (int ell = ellmin; ell != ellmax; ++ell)
        {
            double gdlambda = gaussian(_lambdav[ell], center, sigma) * _dlambdav[ell];
            gsum += gdlambda;
            Jsum += Jv[ell] * gdlambda;
        }
        if (abs(gsum - 1.) > MAX_GAUSS_ERROR_WARN)
        {
            auto units = find<Units>();
            vector<string> message1 = {
                "Integral of Gaussian line profile over radiation field is inaccurate for ",
                " " + _name + " for transition (" + StringUtils::toString(up) + "-" + StringUtils::toString(low) + ")",
                "  integral equals " + StringUtils::toString(gsum) + " rather than unity",
                "  over wavelengths from " + StringUtils::toString(units->owavelength(lambdamin)) + " "
                    + units->uwavelength() + " to " + StringUtils::toString(units->owavelength(lambdamax)) + " "
                    + units->uwavelength() + "."};
            vector<string> message2 = {" 1. Set the wavelength coverage from a velocity window of ±5 x total turbulent "
                                       "velocity (vturb) (i.e., Vmin = -5 vturb, Vmax = +5 vturb) for the radiation "
                                       "field and sample it with around 100 points. The total turbulent velocity "
                                       "includes the micro-turbulent velocity and thermal velocity. Now, vturb = "
                                       + StringUtils::toString(units->ovelocity(sigma)) + " " + units->uvelocity()
                                       + "."};

            if (abs(gsum - 1.) > MAX_GAUSS_ERROR_FAIL && errorForGaussianIntegral())
            {
                log->info("Gausss(" + StringUtils::toString(_lambdav[ellmin])
                          + ")=" + StringUtils::toString(gaussian(_lambdav[ellmin], center, sigma)) + "Gausss("
                          + StringUtils::toString(_lambdav[ellmax - 1])
                          + ")=" + StringUtils::toString(gaussian(_lambdav[ellmax - 1], center, sigma)));
                log->info(message2[0]);
                throw FATALERROR(StringUtils::join(message1, "\n"));
            }

            log->warning(message1[0]);
            log->info(message2[0]);
            for (size_t i = 1; i != message1.size(); ++i) find<Log>()->info(message1[i]);
        }
        // divide by gsum to correct for the fact that, on the discrete wavelength grid, the integral of
        // g over the profile range is only approximately 1 (this is what the check above verifies); if
        // there happen to be no grid points in range at all, gsum and Jsum are both zero and we skip the
        // division to avoid turning that (harmless, J=0) case into a NaN
        double J = gsum > 0. ? Jsum / gsum : Jsum;
        if (storeMeanIntensities()) state->setMeanIntensity(k, J);

        if (!std::isfinite(J))
        {
            throw FATALERROR("Mean intensity J is not finite for transition (" + StringUtils::toString(up) + "-"
                             + StringUtils::toString(low) + ") of " + _name + "The value of J is "
                             + StringUtils::toString(J) + ". The line center is "
                             + StringUtils::toString(_model.center[k]) + ".");
        }

        env.meanJ[k] = J;
    }

    // solve the statistical equilibrium equations with the shared solver; it handles both the
    // radiative terms (Einstein A always, Bul/Blu weighted by env.meanJ) and the collisional terms
    // (using _model.colPartner, with the same robustness against degenerate rates and densities),
    // and throws FatalError directly on a singular matrix or non-finite solution
    vector<double> solution = _gasLineEmission.solveLevelPopulations(_model, env);

    // update the level populations, keeping track of the amount of change
    double change = 0.;
    for (int p = 0; p != numLevels; ++p)
    {
        double oldPop = state->levelPopulation(p);
        double newPop = solution[p];
        state->setLevelPopulation(p, newPop);
        // a population that dropped to exactly zero is a full (100%) relative change; one that
        // stayed at zero is no change -- either way, avoid the 0/0 or x/0 that oldPop/newPop would
        // otherwise produce, which would silently read as "converged" in the caller
        if (newPop > 0.)
            change += abs(oldPop / newPop - 1.);
        else if (oldPop > 0.)
            change += 1.;
    }
    return change / numLevels;
}

////////////////////////////////////////////////////////////////////

bool NonLTELineGasMix::isSpecificStateConverged(int numCells, int /*numUpdated*/, int numNotConverged,
                                                MaterialState* currentAggregate, MaterialState* previousAggregate) const
{
    if (updateDynamicStatesFlag() == false)
    {
        return true;
    }
    else
    {
        // calculate fraction of not converged cells
        double fractionNotConverged = static_cast<double>(numNotConverged) / static_cast<double>(numCells);

        // calculate maximum relative difference between level populations of previous and current iteration
        double changeInGlobalLevelPops = 0.;
        for (int p = 0; p != _model.numLevels(); ++p)
        {
            double currentPop = currentAggregate->levelPopulation(p);
            double previousPop = previousAggregate->levelPopulation(p);
            // as in solveLevelPopulations(), avoid dividing by a previousPop that is exactly zero:
            // treat a population that appeared from zero as a full (100%) change, and one that
            // stayed at zero as no change, instead of letting a NaN silently fail the comparison below
            double diff =
                previousPop > 0. ? abs((currentPop - previousPop) / previousPop) : (currentPop > 0. ? 1. : 0.);
            if (diff > changeInGlobalLevelPops) changeInGlobalLevelPops = diff;
        }

        // log convergence info
        auto log = find<Log>();
        log->info("NonLTELineGasMix convergence info:");
        log->info("  Fraction of not converged cells is " + StringUtils::toString(fractionNotConverged * 100., 'f', 2)
                  + "% (convergence criterion is "
                  + StringUtils::toString(maxFractionNotConvergedCells() * 100., 'f', 2) + "%)");
        log->info("  Global level populations changed by "
                  + StringUtils::toString(changeInGlobalLevelPops * 100., 'f', 2)
                  + "% compared to previous iteration (convergence criterion is "
                  + StringUtils::toString(maxChangeInGlobalLevelPopulations() * 100., 'f', 2) + "%)");

        // convergence is reached when both criteria are satisfied
        return fractionNotConverged <= maxFractionNotConvergedCells()
               && changeInGlobalLevelPops <= maxChangeInGlobalLevelPopulations();
    }
}

////////////////////////////////////////////////////////////////////

double NonLTELineGasMix::mass() const
{
    return _model.mass;
}

////////////////////////////////////////////////////////////////////

double NonLTELineGasMix::sectionAbs(double /*lambda*/) const
{
    return 0.;
}

////////////////////////////////////////////////////////////////////

double NonLTELineGasMix::sectionSca(double /*lambda*/) const
{
    return 0.;
}

////////////////////////////////////////////////////////////////////

double NonLTELineGasMix::sectionExt(double /*lambda*/) const
{
    return 0.;
}

////////////////////////////////////////////////////////////////////

double NonLTELineGasMix::opacityAbs(double lambda, const MaterialState* state, const PhotonPacket* /*pp*/) const
{
    double opacity = 0.;

    // if the cell does not contain any material for this component, leave the opacity at zero
    if (state->numberDensity() > 0.)
    {
        // accumulate the opacities for all radiational transitions
        for (int k = 0; k != _model.numLines(); ++k)
        {
            if (_model.branchRatio[k] < lowestBranchingRatio()) continue;
            double center = _model.center[k];
            double sigma = sigmaForLine(center, state->temperature(), _model.mass);
            Range range(center - PROFILE_RANGE * sigma, center + PROFILE_RANGE * sigma);

            // calculate opacity only if the requested wavelength is in the line profile range
            if (range.contains(lambda))
            {
                int up = _model.indexUpRad[k];
                int low = _model.indexLowRad[k];
                double upnumber = state->levelPopulation(up);
                double lownumber = state->levelPopulation(low);
                double transrate = lownumber * _model.einsteinBlu[k] - upnumber * _model.einsteinBul[k];
                if (transrate != 0. && up < maxUpperLevelForRadiation())
                {
                    constexpr double front = Constants::h() * Constants::c() / 4. / M_PI;
                    opacity += front / center * transrate * gaussian(lambda, center, sigma);
                }
            }
        }

        // apply lower limit to (negative) optical depth
        if (opacity < 0.)
        {
            double diagonal = 1.7320508 * cbrt(state->volume());  // correct only for cubical cell
            if (opacity * diagonal < lowestOpticalDepth()) opacity = lowestOpticalDepth() / diagonal;
        }
    }
    return opacity;
}

////////////////////////////////////////////////////////////////////

double NonLTELineGasMix::opacitySca(double /*lambda*/, const MaterialState* /*state*/, const PhotonPacket* /*pp*/) const
{
    return 0.;
}

////////////////////////////////////////////////////////////////////

double NonLTELineGasMix::opacityExt(double lambda, const MaterialState* state, const PhotonPacket* pp) const
{
    return opacityAbs(lambda, state, pp);
}

////////////////////////////////////////////////////////////////////

bool NonLTELineGasMix::peeloffScattering(double& /*I*/, double& /*Q*/, double& /*U*/, double& /*V*/, double& /*lambda*/,
                                         Direction /*bfkobs*/, Direction /*bfky*/, const MaterialState* /*state*/,
                                         const PhotonPacket* /*pp*/) const
{
    return false;
}

////////////////////////////////////////////////////////////////////

void NonLTELineGasMix::performScattering(double /*lambda*/, const MaterialState* /*state*/, PhotonPacket* /*pp*/) const
{}

////////////////////////////////////////////////////////////////////

Array NonLTELineGasMix::lineEmissionCenters() const
{
    return NR::array(_model.center);
}

////////////////////////////////////////////////////////////////////

Array NonLTELineGasMix::lineEmissionMasses() const
{
    Array masses(_model.numLines());
    for (int k = 0; k != _model.numLines(); ++k) masses[k] = _model.mass;
    return masses;
}

////////////////////////////////////////////////////////////////////

Array NonLTELineGasMix::lineEmissionSpectrum(const MaterialState* state, const Array& /*Jv*/) const
{
    Array luminosities(_model.numLines());
    if (state->numberDensity() > 0.)
    {
        double front = Constants::h() * Constants::c() * state->volume();
        for (int k = 0; k != _model.numLines(); ++k)
        {
            int up = _model.indexUpRad[k];
            if (_model.branchRatio[k] < lowestBranchingRatio() || up >= maxUpperLevelForRadiation())
            {
                luminosities[k] = 0.;
            }
            else
            {
                luminosities[k] = front / _model.center[k] * _model.einsteinA[k] * state->levelPopulation(up);
            }
        }
    }
    return luminosities;
}

////////////////////////////////////////////////////////////////////

double NonLTELineGasMix::indicativeTemperature(const MaterialState* state, const Array& /*Jv*/) const
{
    return state->temperature();
}

////////////////////////////////////////////////////////////////////
