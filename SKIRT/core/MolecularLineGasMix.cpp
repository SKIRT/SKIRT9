/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "MolecularLineGasMix.hpp"
#include "Configuration.hpp"
#include "Constants.hpp"
#include "DisjointWavelengthGrid.hpp"
#include "FatalError.hpp"
#include "Log.hpp"
#include "MaterialState.hpp"
#include "StringUtils.hpp"
#include "TextInFile.hpp"
#include "Units.hpp"

////////////////////////////////////////////////////////////////////

void MolecularLineGasMix::setupSelfBefore()
{
    EmittingGasMix::setupSelfBefore();

    // don't do anything if the simulation has no secondary emission
    auto config = find<Configuration>();
    if (!config->hasSecondaryEmission()) return;

    // get the name of the configured species and the name(s) of its collision partners
    string name;
    vector<string> colNames{"H2"};  // just molecular hydrogen by default
    switch (species())
    {
        case Species::Test: name = "TT"; break;
        case Species::Hydroxyl: name = "OH"; break;
        case Species::Formyl: name = "HCO+"; break;
        case Species::CarbonMonoxide: name = "CO"; break;
        case Species::Carbon:
            name = "C";
            colNames = {"H2", "H", "H+", "e-", "He"};
            break;
    }

    // load the mass of the selected species
    {
        TextInFile infile(this, name + "_Mass.txt", "mass", true);
        infile.addColumn("Molecular weight");
        double weight;
        if (infile.readRow(weight)) _mass = weight * Constants::Mproton();
    }

    // load the energy levels
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

    // calculate the Einstein Bul coefficients and the line centers
    _einsteinBul.resize(_numLines);
    _center.resize(_numLines);
    for (int k = 0; k != _numLines; ++k)
    {
        _center[k] = Constants::h() * Constants::c() / (_energy[_indexUpRad[k]] - _energy[_indexLowRad[k]]);
        _einsteinBul[k] = _einsteinA[k] * pow(_center[k], 5.) / (2. * Constants::h() * Constants::c() * Constants::c());
    }

    // load the collisional transitions for each interaction partner
    bool first = true;
    for (const auto& colName : colNames)
    {
        // add a new empty data structure and get a writable reference to it
        _colPartner.emplace_back();
        auto& partner = _colPartner.back();

        // load the temperature grid points
        {
            TextInFile infile(this, name + "_Col_" + colName + "_Temp.txt", "temperature grid", true);
            infile.addColumn("Temperature", "temperature", "K");
            double temperature;
            while (infile.readRow(temperature)) partner.T.push_back(temperature);
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
    log->info("Radiative lines actually in use for " + name + ":");
    if (_numLines == 0) throw FATALERROR("There are no radiative transitions; increase the number of energy levels");
    for (int k = 0; k != _numLines; ++k)
    {
        log->info("  (" + StringUtils::toString(_indexUpRad[k]) + "-" + StringUtils::toString(_indexLowRad[k]) + ") "
                  + StringUtils::toString(units->owavelength(_center[k])) + " " + units->uwavelength());

        // verify that the radiation field wavelength grid has a bin covering the line center
        if (config->radiationFieldWLG()->bin(_center[k]) < 0)
            throw FATALERROR("Radiation field wavelength grid does not cover the central line for this transition");
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
    return {SnapshotParameter::custom("H2 number density", "numbervolumedensity", "1/cm3"),
            SnapshotParameter::custom("Micro turbulence", "velocity", "km/s")};
}

////////////////////////////////////////////////////////////////////

vector<StateVariable> MolecularLineGasMix::specificStateVariableInfo() const
{
    vector<StateVariable> result{StateVariable::numberDensity(), StateVariable::temperature(),
                                 StateVariable::custom(0, "H2 number density", "numbervolumedensity")};

    // add one custom variable for each level population (the indices start at one because index zero is taken
    // by the H2 number density)
    for (int p = 1; p <= _numEnergyLevels; ++p)
        result.push_back(
            StateVariable::custom(p, "level population " + StringUtils::toString(p), "numbervolumedensity"));
    result.push_back(StateVariable::custom(_numEnergyLevels + 1, "convergence", ""));
    result.push_back(StateVariable::custom(_numEnergyLevels + 2, "Micro turbulence", "velocity"));
    result.push_back(StateVariable::custom(_numEnergyLevels + 3, "temperature", "temperature"));
    for (int p = 1; p <= _numLines; ++p)
        result.push_back(StateVariable::custom(_numEnergyLevels + 3 + p, "radiation field", ""));

    return result;
}

////////////////////////////////////////////////////////////////////

void MolecularLineGasMix::initializeSpecificState(MaterialState* state, double /*metallicity*/, double temperature,
                                                  const Array& params) const
{
    (void)state;
    (void)temperature;
    (void)params;
    /*
    // leave the properties untouched if the cell does not contain any material for this component
    if (state->numberDensity() > 0.)
    {
        // collisional partner density
        state->setCustom(0, params.size() ? params[0] : state->numberDensity() * defaultCollisionPartnerRatios()[0]);

        // kinetic temperature
        double Tkin = temperature >= 0. ? temperature : defaultTemperature();
        state->setCustom(numLevelPops + 3, Tkin);

        // effective temperature (including micro-turbulence)
        double vturb = params.size() ? params[1] : defaultMicroTurbulenceVelocity();
        double Teff = Tkin + vturb * vturb * mass() / Constants::k();
        state->setTemperature(Teff);

        // initialization of the level population using boltzmann distribution (start with LTE)
        double normN = 0.;
        Array boltzDist(numLevelPops);
        for (int p = 0; p < numLevelPops; ++p)
        {
            boltzDist[p] = _weights[p] * exp(-_energyStates[p] / Constants::k() / state->custom(numLevelPops + 3));
            normN += boltzDist[p];
        }
        for (int p = 0; p < numLevelPops; ++p)
        {
            double statePop = state->numberDensity() * boltzDist[p] / normN;
            state->setCustom(p + 1, statePop);
        }
    }
*/
}

////////////////////////////////////////////////////////////////////
/*
namespace
{
    // return a collisional de-excitation coefficeint using log-linear interpolation based on the table of
    // the collisional coefficients
    double coefficientCul(double temperature, int indexTrans, int indexTem, int numTemp, const Array& tableTemp,
                          const vector<vector<double>>& tableCul)
    {
        double Cul = 0.;
        if (indexTem < 0)
        {
            throw FATALERROR("Failed calculating a Collisional coefficient of transition "
                             + StringUtils::toString(indexTrans));
        }
        else if (indexTem == 0)
            Cul = tableCul[indexTrans][0];
        else if (indexTem == numTemp - 1)
            Cul = tableCul[indexTrans][numTemp - 1];
        else
        {
            double logTemp = log(temperature);
            double logTempTableHigh = log(tableTemp[indexTem + 1]);
            double logTempTableLow = log(tableTemp[indexTem]);
            double logCulHigh = log(tableCul[indexTrans][indexTem + 1]);
            double logCulLow = log(tableCul[indexTrans][indexTem]);
            double steepness = (logCulHigh - logCulLow) / (logTempTableHigh - logTempTableLow);
            Cul = exp(steepness * (logTemp - logTempTableLow) + logCulLow);
        }
        return Cul;
    }

    // return a collisional excitation coefficeint Clu using a collisional de-excitation coefficeint Cul
    double coefficientClu(double Cul, double temperature, int up, int low, const vector<double>& weight,
                          const vector<double>& energy)
    {
        double weightRatio = weight[up] / weight[low];
        double dif_energy = abs(energy[up] - energy[low]);
        double expTerm = 0.;

        expTerm = exp(-dif_energy / Constants::k() / temperature);
        double Clu = Cul * weightRatio * expTerm;
        return Clu;
    }

    // solve the given inverse matrix using LU decomposition
    Array solveInverseMatrix(vector<vector<double>> matrixPop)
    {
        size_t rowSize = matrixPop.size();
        size_t collumnSize = matrixPop[0].size();
        Array solution(rowSize);
        double inverseMat;

        // forwarding elimination
        for (size_t i = 0; i < rowSize; i++) solution[i] = matrixPop[i][collumnSize - 1];

        for (size_t k = 0; k < rowSize - 1; ++k)
        {
            if (matrixPop[k][k] == 0.0)
            {
                std::swap(matrixPop[k], matrixPop[k + 1]);
                std::swap(solution[k], solution[k + 1]);
            }
            for (size_t i = k + 1; i < collumnSize - 1; ++i)
            {
                inverseMat = matrixPop[i][k] / matrixPop[k][k];
                for (size_t j = k + 1; j < collumnSize - 1; ++j) matrixPop[i][j] -= inverseMat * matrixPop[k][j];
                matrixPop[i][k] = inverseMat;
            }
        }

        // forwarding elimination
        for (size_t i = 0; i < rowSize; ++i)
            for (size_t j = 0; j < i; ++j) solution[i] -= matrixPop[i][j] * solution[j];

        // backward substitution
        for (int i = static_cast<int>(rowSize) - 1; i >= 0; --i)
        {
            for (size_t j = i + 1; j < collumnSize - 1; ++j) solution[i] -= matrixPop[i][j] * solution[j];
            solution[i] /= matrixPop[i][i];
        }

        // return the solution of given inverse matrix
        return solution;
    }

    double m_SQRT_2PI = M_SQRT1_2 * sqrt(M_1_PI);

    // return gaussian profile using 4 parameters: wavelength, wavelength of the line center,
    // temperature, molecular mass
    double gaussianProfile(double lambda, double center, double const_sigmatherm)
    {
        double sigmatherm = center * const_sigmatherm;
        double front = 1. / (sigmatherm)*m_SQRT_2PI;
        double back = exp(-(lambda - center) * (lambda - center) / 2.0 / sigmatherm / sigmatherm);
        return front * back;
    }

    // this number divide the dispersion of gaussian (sigma) to determin differential wavelength for the intergration
    //of the gaussian profile
    double divnum = 5.0;
    // this number determin the range of wavelength for the intergration of the gaussian profile
    int rangeLambdaNum = 20;  //if:=25 (5sigma);error 1e-4%, if:=20 (4sigma); error: 1e-2%, if:=15 (3sigma); error: 0.3%

}
*/
////////////////////////////////////////////////////////////////////

UpdateStatus MolecularLineGasMix::updateSpecificState(MaterialState* state, const Array& Jv) const
{
    UpdateStatus status;
    (void)state;
    (void)Jv;
    /*
    // leave the properties untouched if the cell does not contain any material for this component
    if (state->numberDensity() > 0.0)
    {
        auto config = find<Configuration>();

        double numHydroMol = state->custom(0);

        Array levelpop_before(numLevelPops);
        for (int p = 1; p <= numLevelPops; ++p) levelpop_before[p - 1] = state->custom(p);

        // set up mean Intensity
        Array rf(_numLines);

        // wavelength of transition lines
        Array centers = lineEmissionCenters();
        double const_sigmatherm = 1.0 / Constants::c() * sqrt(Constants::k() * state->temperature() / mass());

        // calculate radiation fields at each transition line using the line profile
        for (int i = 0; i < _numLines; ++i)
        {

            // integrate radiation J fields with a weight of gaussian profile
            double difLambda = centers[i] * const_sigmatherm / divnum;
            double lambdaNow = centers[i] - difLambda * rangeLambdaNum;
            double lambdaEnd = centers[i] + difLambda * rangeLambdaNum;

            rf[i] = 0.;

            double integralsum = 0.;
            int count_i = 0;
            while (lambdaNow < lambdaEnd)
            {
                double lambdaPre = lambdaNow;
                int indexJv = -1;
                indexJv = config->radiationFieldWLG()->bin(lambdaNow);
                if (indexJv >= 0 and indexJv < config->radiationFieldWLG()->numBins())
                {
                    double rightJvBorder = config->radiationFieldWLG()->rightBorder(indexJv);
                    lambdaNow += difLambda;
                    if (lambdaNow > rightJvBorder)
                    {
                        double gaussValue =
                            gaussianProfile((rightJvBorder + lambdaPre) / 2., centers[i], const_sigmatherm)
                            * (rightJvBorder - lambdaPre);
                        rf[i] += Jv[indexJv] * gaussValue;
                        integralsum += gaussValue;
                        gaussValue = gaussianProfile((rightJvBorder + lambdaNow) / 2., centers[i], const_sigmatherm)
                                     * (lambdaNow - rightJvBorder);
                        if (indexJv + 1 < config->radiationFieldWLG()->numBins())
                        {
                            rf[i] += Jv[indexJv + 1] * gaussValue;
                            integralsum += gaussValue;
                        }
                        else
                            break;
                    }
                    else
                    {
                        double gaussValue =
                            difLambda * gaussianProfile((lambdaNow + lambdaPre) / 2., centers[i], const_sigmatherm);
                        rf[i] += Jv[indexJv] * gaussValue;
                        integralsum += gaussValue;
                    }
                }
                else if (indexJv < 0)
                    lambdaNow += difLambda;
                else
                {
                    find<Log>()->info("Break Loop in calculation of radiation fields, index "
                                      + StringUtils::toString(indexJv) + " and count "
                                      + StringUtils::toString(count_i));
                    break;
                }
                count_i++;
            }

            if (integralsum > 1.01 or integralsum < 0.99)
                find<Log>()->warning(
                    "Integrated Gaussian = " + StringUtils::toString(integralsum)
                    + ", please adjust lineRange in updateSpecificState, wavelengthgrids, or temperature");
        }

        // statistical equilibrium matrix for level population
        vector<vector<double>> matrixPop(numLevelPops, vector<double>(numLevelPops + 1, 0.));

        // temperature of the carbon monoxide to determine level populations
        double tem = state->custom(numLevelPops + 3);

        // the most closet index of the templeature in the collisional table to the temperature of the state
        int indexTem = 0;

        // add the terms of the raditational transtions
        for (int i = 0; i < _numLines; ++i)
        {
            int up = _indexUpRad[_indexRadTrans[i]];
            int low = _indexLowRad[_indexRadTrans[i]];

            // add Eistein Aul coefficeints (spontaneous emission) in the statistical equibrium matrix
            matrixPop[up][up] += _einsteinA[_indexRadTrans[i]];
            matrixPop[low][up] -= _einsteinA[_indexRadTrans[i]];

            // add Eistein Bul coefficeints (stimulated emission) in the statistical equibrium matrix
            matrixPop[up][up] += _einsteinBul[_indexRadTrans[i]] * rf[i];
            matrixPop[low][up] -= _einsteinBul[_indexRadTrans[i]] * rf[i];

            // add Eistein Blu coefficeints (absorption) in the statistical equibrium matrix
            double wRatio = _weights[up] / _weights[low];
            matrixPop[low][low] += wRatio * _einsteinBul[_indexRadTrans[i]] * rf[i];
            matrixPop[up][low] -= wRatio * _einsteinBul[_indexRadTrans[i]] * rf[i];
        }

        // add the terms of the collisional transtions with hydorgen molecules
        for (size_t i = 0; i < _indexColTrans.size(); ++i)
        {
            int up = _indexUpCol[_indexColTrans[i]];
            int low = _indexLowCol[_indexColTrans[i]];
            double Cul =
                coefficientCul(tem, _indexColTrans[i], indexTem, numTempereature, temperatureTable, _collisionCul);
            double Clu = coefficientClu(Cul, tem, up, low, _weights, _energyStates);

            matrixPop[up][up] += Cul * numHydroMol;
            matrixPop[low][low] += Clu * numHydroMol;
            matrixPop[up][low] -= Clu * numHydroMol;
            matrixPop[low][up] -= Cul * numHydroMol;
        }

        // add the normalization of number density in the last row of the matrix.
        for (int j = 0; j <= numLevelPops; ++j)
        {
            if (j == numLevelPops)
                matrixPop[numLevelPops - 1][j] = state->numberDensity();
            else
                matrixPop[numLevelPops - 1][j] = 1.;
        }

        // solve the given inverse matrix
        Array solution = solveInverseMatrix(matrixPop);

        // update the level population
        double total_num = 0.;
        for (int p = 1; p <= numLevelPops; ++p)
        {
            if (solution[p - 1] >= state->numberDensity())
                throw FATALERROR(StringUtils::toString(p - 1)
                                 + " level: Failed updating level population, the number density in the level ="
                                 + StringUtils::toString(solution[p - 1]));
            // set level population
            else
                state->setCustom(p, solution[p - 1]);
            total_num += solution[p - 1];
        }

        // for debug: confirm number density
        if (total_num >= state->numberDensity() * (1.0 + 1e-5) || total_num <= state->numberDensity() * (1.0 - 1e-5))
            throw FATALERROR("Failed Update total: " + StringUtils::toString(state->numberDensity())
                             + "; after: " + StringUtils::toString(total_num));
        // confirm the convergence of level population
        double converge_value = 0.;
        for (int p = 1; p <= numLevelPops; ++p)
            if (p < numLevelPops) converge_value += abs((levelpop_before[p - 1] / state->custom(p) - 1));
        converge_value /= (numLevelPops - 1);
        if (converge_value > maxChangeInLevelPopulations())
            status.updateNotConverged();
        else
            status.updateConverged();

        // memorize the radiation fields at line centers
        for (int p = 1; p <= _numLines; ++p) state->setCustom(numLevelPops + 3 + p, rf[p - 1]);
    }
*/
    return status;
}

////////////////////////////////////////////////////////////////////

bool MolecularLineGasMix::isSpecificStateConverged(int numCells, int /*numUpdated*/, int numNotConverged) const
{
    return static_cast<double>(numNotConverged) / static_cast<double>(numCells) <= maxFractionUnconvergedCells();
}

////////////////////////////////////////////////////////////////////

double MolecularLineGasMix::mass() const
{
    return _mass;
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
    double kappa = 0.;
    (void)lambda;
    (void)state;
    /*
    if (state->numberDensity() > 0.0)
    {
        constexpr double h = Constants::h();
        constexpr double c = Constants::c();
        Array numPop(numLevelPops);
        for (int i = 0; i < numLevelPops; ++i) numPop[i] = state->custom(i + 1);

        // wavelength range outside of which we calculate, the line opacities are zero (range of plus-minus lambda)
        Range lineRange(0.97 * lambda, 1.03 * lambda);

        // wavelengths of the rotational lines
        Array centers = lineEmissionCenters();
        double const_sigmatherm = 1.0 / Constants::c() * sqrt(Constants::k() * state->temperature() / mass());

        // calculate the total opacity including all of the rotational transition lines
        for (int p = 0; p < _numLines; ++p)
        {
            if (lineRange.contains(centers[p]))
            {
                int upindex = _indexUpRad[_indexRadTrans[p]];
                int lowindex = _indexLowRad[_indexRadTrans[p]];
                double upnumber = numPop[upindex];
                double lownumber = numPop[lowindex];
                if (upnumber != 0.0 && lownumber != 0.0)
                {
                    double einsteinBlu = _weights[upindex] / _weights[lowindex] * _einsteinBul[_indexRadTrans[p]];
                    double transrate = lownumber * einsteinBlu - upnumber * _einsteinBul[_indexRadTrans[p]];
                    kappa += h * c / 4.0 / M_PI / centers[p] * transrate
                             * gaussianProfile(lambda, centers[p], const_sigmatherm);
                }
            }
        }
    }
*/
    return kappa;
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
    return _center;
}

////////////////////////////////////////////////////////////////////

Array MolecularLineGasMix::lineEmissionMasses() const
{
    Array masses(_numLines);
    for (int p = 0; p < _numLines; ++p) masses[p] = mass();
    return masses;
}

////////////////////////////////////////////////////////////////////

Array MolecularLineGasMix::lineEmissionSpectrum(const MaterialState* state, const Array& /*Jv*/) const
{
    Array luminosities(_numLines);
    (void)state;
    /*
    if (state->numberDensity() > 0.0)
    {
        Array centers = lineEmissionCenters();
        Array numPop(numLevelPops);
        for (int i = 0; i < numLevelPops; ++i) numPop[i] = state->custom(i + 1);

        //caluculate the line luminosities in the cell
        for (int p = 0; p < _numLines; ++p)
        {
            luminosities[p] = Constants::h() * Constants::c() / centers[p] * _einsteinA[_indexRadTrans[p]]
                              * state->custom(_indexUpRad[_indexRadTrans[p]] + 1) * state->volume();
        }
    }
*/
    return luminosities;
}

////////////////////////////////////////////////////////////////////

double MolecularLineGasMix::indicativeTemperature(const MaterialState* state, const Array& /*Jv*/) const
{
    return state->temperature();
}

////////////////////////////////////////////////////////////////////
