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
        case Species::Test: name = "TT"; break;
        case Species::Hydroxyl: name = "OH"; break;
        case Species::Formyl: name = "HCO+"; break;
        case Species::CarbonMonoxide: name = "CO"; break;
        case Species::AtomicCarbon:
            name = "C";
            colNames = {"H2", "H", "H+", "e-", "He"};
            break;
        case Species::IonizedCarbon:
            name = "C+";
            colNames = {"H2", "H", "e-"};
            break;
    }

    // load the mass of the selected species
    {
        TextInFile infile(this, name + "_Mass.txt", "mass", true);
        infile.addColumn("Mass", "mass", "amu");
        double weight;
        if (infile.readRow(weight)) _mass = weight;
    }

    // load the energy levels and weights
    {
        TextInFile infile(this, name + "_Energy.txt", "energy levels", true);
        infile.addColumn("Energy", "energy", "1/cm");
        infile.addColumn("Weight");
        double energy, weight;
        while (infile.readRow(energy, weight))
        {
            _energy.push_back(energy);
            _weight.push_back(weight);
            if (_energy.size() == static_cast<size_t>(numEnergyLevels())) break;
        }
        _numLevels = _energy.size();
    }

    // load the radiative transitions
    {
        TextInFile infile(this, name + "_Rad_Coeff.txt", "radiative transitions", true);
        infile.addColumn("Up index");
        infile.addColumn("Low index");
        infile.addColumn("Einstein A", "transitionrate", "1/s");
        double up, low, rate;
        while (infile.readRow(up, low, rate))
        {
            int indexUp = up;
            int indexLow = low;
            if (indexUp < _numLevels && indexLow < _numLevels)
            {
                _indexUpRad.push_back(indexUp);
                _indexLowRad.push_back(indexLow);
                _einsteinA.push_back(rate);
            }
        }
    }
    _numLines = _indexUpRad.size();

    // calculate the line centers and the Einstein B coefficients
    _center.resize(_numLines);
    _einsteinBul.resize(_numLines);
    _einsteinBlu.resize(_numLines);
    for (int k = 0; k != _numLines; ++k)
    {
        _center[k] = Constants::h() * Constants::c() / (_energy[_indexUpRad[k]] - _energy[_indexLowRad[k]]);
        _einsteinBul[k] = _einsteinA[k] * pow(_center[k], 5.) / (2. * Constants::h() * Constants::c() * Constants::c());
        _einsteinBlu[k] = _einsteinBul[k] * _weight[_indexUpRad[k]] / _weight[_indexLowRad[k]];
    }

    // load the collisional transitions for each interaction partner
    bool first = true;
    for (const auto& colName : colNames)
    {
        // add a new empty data structure and get a writable reference to it
        _colPartner.emplace_back();
        auto& partner = _colPartner.back();
        partner.name = colName;

        // load the temperature grid points
        {
            TextInFile infile(this, name + "_Col_" + colName + "_Temp.txt", "temperature grid", true);
            infile.addColumn("Temperature", "temperature", "K");
            infile.readAllColumns(partner.T);
        }

        // load the transition indices (just for the first partner) and coefficients (for every partner)
        {
            int numTemperatures = partner.T.size();
            TextInFile infile(this, name + "_Col_" + colName + "_Coeff.txt", "collisional transitions", true);
            infile.addColumn("Up index");
            infile.addColumn("Low index");
            for (int i = 0; i != numTemperatures; ++i) infile.addColumn("Collisional K", "collisionalrate", "cm3/s");
            Array row;
            while (infile.readRow(row))
            {
                int indexUp = row[0];
                int indexLow = row[1];
                if (indexUp < _numLevels && indexLow < _numLevels)
                {
                    if (first)
                    {
                        _indexUpCol.push_back(indexUp);
                        _indexLowCol.push_back(indexLow);
                    }
                    Array coeff(numTemperatures);
                    for (int i = 0; i != numTemperatures; ++i) coeff[i] = row[i + 2];
                    partner.Kul.emplace_back(coeff);
                }
            }
        }
        first = false;
    }
    _numColTrans = _indexUpCol.size();
    _numColPartners = colNames.size();

    // log summary info on the radiative lines
    auto log = find<Log>();
    auto units = find<Units>();
    log->info("Radiative lines for " + name + ":");
    if (_numLines == 0) throw FATALERROR("There are no radiative transitions; increase the number of energy levels");
    for (int k = 0; k != _numLines; ++k)
    {
        log->info("  (" + StringUtils::toString(_indexUpRad[k]) + "-" + StringUtils::toString(_indexLowRad[k]) + ") "
                  + StringUtils::toString(units->owavelength(_center[k])) + " " + units->uwavelength());
    }

    // log summary info on the collisional partner(s)
    log->info("Collisional partner(s) for " + name + ": " + StringUtils::join(colNames, ", "));

    // verify that the radiation field wavelength grid, if present, has a bin covering the line centers
    // and cache the characteristic wavelengths and bin widths
    auto rfwlg = find<Configuration>()->radiationFieldWLG();
    if (rfwlg)
    {
        rfwlg->setup();
        for (int k = 0; k != _numLines; ++k)
        {
            if (rfwlg->bin(_center[k]) < 0)
                throw FATALERROR("Radiation field wavelength grid does not cover the central line for transition ("
                                 + StringUtils::toString(_indexUpRad[k]) + "-" + StringUtils::toString(_indexLowRad[k])
                                 + ")");
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
        for (int p = 0; p != _numLevels; ++p)
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
    for (const auto& partner : _colPartner)
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
    for (const auto& partner : _colPartner)
        result.push_back(StateVariable::custom(index++, partner.name + " number density", "numbervolumedensity"));

    // add custom variable for the population of each energy level
    const_cast<NonLTELineGasMix*>(this)->_indexFirstLevelPopulation = index;
    for (int p = 0; p != _numLevels; ++p)
        result.push_back(
            StateVariable::custom(index++, "population of level " + std::to_string(p), "numbervolumedensity"));

    // if requested, add custom variable for the line-profile-averaged mean intensity at each transition line
    if (storeMeanIntensities())
    {
        const_cast<NonLTELineGasMix*>(this)->_indexFirstMeanIntensity = index;
        for (int k = 0; k != _numLines; ++k)
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
        // copy kinetic temperature from import or default
        double Tkin = temperature >= 0. ? temperature : defaultTemperature();
        state->setKineticTemperature(Tkin);

        // set effective temperature, including imported or default turbulence
        double vturb = params.size() ? params[_numColPartners] : defaultTurbulenceVelocity();
        double Teff = Tkin + 0.5 * vturb * vturb * _mass / Constants::k();
        state->setTemperature(Teff);

        // copy collisional partner densities from import or default
        if (params.size())
        {
            for (int c = 0; c != _numColPartners; ++c) state->setColPartnerDensity(c, params[c]);
        }
        else
        {
            const auto& ratios = defaultCollisionPartnerRatios();
            if (ratios.size() < _colPartner.size())
                throw FATALERROR("The number of collision partners exceeds the number of default ratios");
            for (int c = 0; c != _numColPartners; ++c)
                state->setColPartnerDensity(c, state->numberDensity() * ratios[c]);
        }

        // initialize level population using boltzmann distribution (i.e., start with LTE)
        Array levelPops(_numLevels);
        for (int p = 0; p != _numLevels; ++p) levelPops[p] = _weight[p] * exp(-_energy[p] / Constants::k() / Tkin);

        // if the user configured a file with initial level populations, use those data instead
        size_t m = state->cellIndex();
        if (m < _initLevelPops.size())
            for (int p = 0; p != _numLevels; ++p) levelPops[p] = _initLevelPops[m][p + 1];

        // normalize and store
        levelPops *= state->numberDensity() / levelPops.sum();
        for (int p = 0; p != _numLevels; ++p) state->setLevelPopulation(p, levelPops[p]);
    }
}

////////////////////////////////////////////////////////////////////

namespace
{
    // solve the square set of linear equations represented by the given matrix using LU decomposition
    // the matrix should have N rows and N+1 columns; its contents is overwritten and
    // the solution is returned as an array of size N
    Array solveMatrixEquation(vector<vector<double>>& matrix)
    {
        size_t size = matrix.size();
        Array solution(size);

        // forwarding elimination
        for (size_t i = 0; i < size; i++) solution[i] = matrix[i][size];

        // decomposition
        for (size_t k = 0; k < size - 1; ++k)
        {
            if (matrix[k][k] == 0.0)
            {
                std::swap(matrix[k], matrix[k + 1]);
                std::swap(solution[k], solution[k + 1]);
            }
            for (size_t i = k + 1; i < size; ++i)
            {
                double inverse = matrix[i][k] / matrix[k][k];
                for (size_t j = k + 1; j < size; ++j) matrix[i][j] -= inverse * matrix[k][j];
                matrix[i][k] = inverse;
            }
        }

        // forwarding elimination
        for (size_t i = 0; i < size; ++i)
            for (size_t j = 0; j < i; ++j) solution[i] -= matrix[i][j] * solution[j];

        // backward substitution
        for (int i = static_cast<int>(size) - 1; i >= 0; --i)
        {
            for (size_t j = i + 1; j < size; ++j) solution[i] -= matrix[i][j] * solution[j];
            solution[i] /= matrix[i][i];
        }

        // return the solution
        return solution;
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
    if (state->numberDensity() > 0)
    {
        // allocate the statistical equilibrium matrix for the level populations
        vector<vector<double>> matrix(_numLevels, vector<double>(_numLevels + 1));

        // add the terms for the radiational transitions
        for (int k = 0; k != _numLines; ++k)
        {
            int up = _indexUpRad[k];
            int low = _indexLowRad[k];

            // add the Einstein Aul coefficients (spontaneous emission)
            matrix[up][up] -= _einsteinA[k];
            matrix[low][up] += _einsteinA[k];

            // calculate the mean intensity of the radiation field convolved over the normalized line profile g:
            //   J_convolved = \int J_lambda(lambda) g(lambda) d lambda  /  \int g(lambda) d lambda
            // we use all wavelength points within a given range around the line center and verify that the
            // grid is sufficiently resolved to reproduce the normalizaton value of 1 = \int g(lambda) d lambda
            double center = _center[k];
            double sigma = sigmaForLine(center, state->temperature(), _mass);
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
                vector<string> message = {
                    "Integral of Gaussian line profile over radiation field is inaccurate:",
                    "  integral equals " + StringUtils::toString(gsum) + " rather than unity",
                    "  over wavelengths from " + StringUtils::toString(units->owavelength(lambdamin)) + " "
                        + units->uwavelength() + " to " + StringUtils::toString(units->owavelength(lambdamax)) + " "
                        + units->uwavelength(),
                    "  --> increase the resolution of the radiation field wavelength grid"};
                if (abs(gsum - 1.) > MAX_GAUSS_ERROR_FAIL) throw FATALERROR(StringUtils::join(message, "\n"));
                auto log = find<Log>();
                log->warning(message[0]);
                for (size_t i = 1; i != message.size(); ++i) find<Log>()->info(message[i]);
            }
            double J = Jsum / gsum;
            if (storeMeanIntensities()) state->setMeanIntensity(k, J);

            // add the Einstein Bul coefficients (stimulated emission)
            matrix[up][up] -= _einsteinBul[k] * J;
            matrix[low][up] += _einsteinBul[k] * J;

            // add the Einstein Blu coefficients (absorption)
            matrix[low][low] -= _einsteinBlu[k] * J;
            matrix[up][low] += _einsteinBlu[k] * J;
        }

        // add the terms for the collisional transitions
        double T = state->kineticTemperature();
        for (int t = 0; t != _numColTrans; ++t)
        {
            int up = _indexUpCol[t];
            int low = _indexLowCol[t];
            double weightRatio = _weight[up] / _weight[low];
            double energyDiff = _energy[up] - _energy[low];
            double Kconversion = weightRatio * exp(-energyDiff / Constants::k() / T);

            // loop over the collisional partners
            for (int c = 0; c != _numColPartners; ++c)
            {
                // determine Kul by interpolation from the temperature-dependent table
                double Kul = NR::clampedValue<NR::interpolateLogLog>(T, _colPartner[c].T, _colPartner[c].Kul[t]);

                // determine Klu from Kul
                double Klu = Kul * Kconversion;

                // add the coefficients after multiplication by the partner number density
                double n = state->colPartnerDensity(c);
                matrix[up][up] -= Kul * n;
                matrix[low][low] -= Klu * n;
                matrix[up][low] += Klu * n;
                matrix[low][up] += Kul * n;
            }
        }

        // replace the last row of the matrix by the normalization of the number density
        for (int p = 0; p != _numLevels; ++p) matrix[_numLevels - 1][p] = 1.;
        matrix[_numLevels - 1][_numLevels] = state->numberDensity();

        // solve the set of equations represented by the matrix
        Array solution = solveMatrixEquation(matrix);

        // update the level populations, keeping track of the amount of change
        double change = 0.;
        for (int p = 0; p != _numLevels; ++p)
        {
            double oldPop = state->levelPopulation(p);
            double newPop = solution[p];
            state->setLevelPopulation(p, newPop);
            change += abs(oldPop / newPop - 1.);
        }
        change /= _numLevels;

        // verify convergence
        if (change > maxChangeInLevelPopulations())
            status.updateNotConverged();
        else
            status.updateConverged();
    }
    return status;
}

////////////////////////////////////////////////////////////////////

bool NonLTELineGasMix::isSpecificStateConverged(int numCells, int /*numUpdated*/, int numNotConverged) const
{
    return static_cast<double>(numNotConverged) / static_cast<double>(numCells) <= maxFractionNotConvergedCells();
}

////////////////////////////////////////////////////////////////////

double NonLTELineGasMix::mass() const
{
    return _mass;
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
        for (int k = 0; k != _numLines; ++k)
        {
            double center = _center[k];
            double sigma = sigmaForLine(center, state->temperature(), _mass);
            Range range(center - PROFILE_RANGE * sigma, center + PROFILE_RANGE * sigma);

            // calculate opacity only if the requested wavelength is in the line profile range
            if (range.contains(lambda))
            {
                int up = _indexUpRad[k];
                int low = _indexLowRad[k];
                double upnumber = state->levelPopulation(up);
                double lownumber = state->levelPopulation(low);
                double transrate = lownumber * _einsteinBlu[k] - upnumber * _einsteinBul[k];
                if (transrate != 0.)
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

void NonLTELineGasMix::peeloffScattering(double& /*I*/, double& /*Q*/, double& /*U*/, double& /*V*/, double& /*lambda*/,
                                         Direction /*bfkobs*/, Direction /*bfky*/, const MaterialState* /*state*/,
                                         const PhotonPacket* /*pp*/) const
{}

////////////////////////////////////////////////////////////////////

void NonLTELineGasMix::performScattering(double /*lambda*/, const MaterialState* /*state*/, PhotonPacket* /*pp*/) const
{}

////////////////////////////////////////////////////////////////////

Array NonLTELineGasMix::lineEmissionCenters() const
{
    return _center;
}

////////////////////////////////////////////////////////////////////

Array NonLTELineGasMix::lineEmissionMasses() const
{
    Array masses(_numLines);
    for (int k = 0; k != _numLines; ++k) masses[k] = _mass;
    return masses;
}

////////////////////////////////////////////////////////////////////

Array NonLTELineGasMix::lineEmissionSpectrum(const MaterialState* state, const Array& /*Jv*/) const
{
    Array luminosities(_numLines);
    if (state->numberDensity() > 0.)
    {
        double front = Constants::h() * Constants::c() * state->volume();
        for (int k = 0; k != _numLines; ++k)
            luminosities[k] = front / _center[k] * _einsteinA[k] * state->levelPopulation(_indexUpRad[k]);
    }
    return luminosities;
}

////////////////////////////////////////////////////////////////////

double NonLTELineGasMix::indicativeTemperature(const MaterialState* state, const Array& /*Jv*/) const
{
    return state->temperature();
}

////////////////////////////////////////////////////////////////////
