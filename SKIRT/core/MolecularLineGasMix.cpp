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

////////////////////////////////////////////////////////////////////

namespace
{
    // number of level populations which you use in this calculation. Users can define this number.
    constexpr int numLevelPops = 2;

    // molecular weight of test molecule for a benchmark 1 of van Zadelho et al. (2002)
    constexpr double molecularWeight = 7.337241790008474;

    // name of the file which contatins energies and weights of energy states
    string filenameEnergy("Test_Energy_States.txt");

    // number of energy levels in the table
    constexpr int numLevelTable = 2;

    // name of the file which contatins up and low indexes and Einstein A coefficients of radiational transitions
    string filenameRadiative("Test_Radiative_Trans.txt");

    // number of radiational transition lines in the table
    constexpr int numLinesTable = 1;

    // name of the file which contatins up and low indexes and Collisional excitation coefficients of
    // collisional transitions
    string filenameCollisional("Test_Collisional_Trans.txt");

    // number of coliisional transition lines in the table
    constexpr int numCollisionTable = 1;

    // temeratures in the table of C coefficients
    const Array temperatureTable = {20.0};

    // number of temperatures included in the table of C coefficients
    constexpr int numTempereature = 1;
}

////////////////////////////////////////////////////////////////////

void MolecularLineGasMix::setupSelfBefore()
{
    EmittingGasMix::setupSelfBefore();

    auto config = find<Configuration>();
    if (config->hasSecondaryEmission())
    {
        // import the energies and the weight of energy states from a text file
        TextInFile infile1(this, filenameEnergy, "energy states", true);
        infile1.addColumn("Energy", "energy", "1/cm");
        infile1.addColumn("Weight");
        auto result1 = infile1.readAllColumns();
        for (int i = 0; i < numLevelTable; ++i)
        {
            _energyStates.push_back(result1[0][i]);
            _weights.push_back(result1[1][i]);
        }
        infile1.close();

        // import the transition indices and Einstein A coefficients of radiational transitions from a text file
        TextInFile infile2(this, filenameRadiative, "radiative transitions", true);
        infile2.addColumn("Up index");
        infile2.addColumn("Low index");
        infile2.addColumn("Einstein A", "transitionrate", "1/s");
        auto result2 = infile2.readAllColumns();
        for (int i = 0; i < numLinesTable; ++i)
        {
            _indexUpRad.push_back(static_cast<int>(result2[0][i]));
            _indexLowRad.push_back(static_cast<int>(result2[1][i]));
            _einsteinA.push_back(result2[2][i]);
        }
        infile2.close();

        // import the transition indices and collisional K coefficients of collisional transitions from a text file
        TextInFile infile3(this, filenameCollisional, "collisional transitions", true);
        infile3.addColumn("Up index");
        infile3.addColumn("Low index");
        for (int i = 0; i < numTempereature; ++i) infile3.addColumn("Collisional K", "collisionalrate", "cm3/s");
        auto result3 = infile3.readAllColumns();
        for (int i = 0; i < numCollisionTable; ++i)
        {
            _indexUpCol.push_back(static_cast<int>(result3[0][i]));
            _indexLowCol.push_back(static_cast<int>(result3[1][i]));
            _collisionCul.push_back(vector<double>());
            for (int j = 0; j < numTempereature; ++j) _collisionCul[i].push_back(result3[2 + j][i]);
        }
        infile3.close();

        // calculate Einstein Bul using Einstein A
        Array centers(numLinesTable);  // wavelengths of the line center of the rotational transitions
        for (int i = 0; i < numLinesTable; ++i)
        {
            centers[i] =
                Constants::h() * Constants::c() / (_energyStates[_indexUpRad[i]] - _energyStates[_indexLowRad[i]]);
            _einsteinBul.push_back(pow(centers[i], 5.0) / 2.0 / Constants::h() / Constants::c() / Constants::c()
                                   * _einsteinA[i]);
        }

        // store indexes of the radiational transitions which you use in this calculation
        for (int i = 0; i < numLinesTable; ++i)
        {
            if (_indexUpRad[i] < numLevelPops && _indexLowRad[i] < numLevelPops) _indexRadTrans.push_back(i);
        }

        for (int i = 0; i < numCollisionTable; ++i)
        {
            if (_indexUpCol[i] < numLevelPops && _indexLowCol[i] < numLevelPops) _indexColTrans.push_back(i);
        }

        // number of the radiational transition lines which you use in this calculation
        _numLines = static_cast<int>(_indexRadTrans.size());
        if (_numLines < 0)
            throw FATALERROR("There is no rotational lines, please increase the number of the energy states");

        // remember the wavelength indexes of the radiation field bin corresponding to the rotational lines.
        auto log = find<Log>();
        for (int i = 0; i < _numLines; ++i)
        {
            log->info("Test molecule (" + StringUtils::toString(_indexUpRad[_indexRadTrans[i]]) + "-"
                      + StringUtils::toString(_indexLowRad[_indexRadTrans[i]]) + ") "
                      + StringUtils::toString(centers[_indexRadTrans[i]] * 1e6) + " [um]");
            _indexLines.push_back(-1);  // initialization of the indexes of the rotational lines
            _indexLines[i] = config->radiationFieldWLG()->bin(centers[_indexRadTrans[i]]);
            if (_indexLines[i] < 0)
                throw FATALERROR("Radiation field wavelength grid does not include a transition "
                                 + StringUtils::toString(_indexRadTrans[i]) + " line");
        }
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
    for (int p = 1; p <= numLevelPops; ++p)
        result.push_back(
            StateVariable::custom(p, "level population " + StringUtils::toString(p), "numbervolumedensity"));
    result.push_back(StateVariable::custom(numLevelPops + 1, "convergence", ""));
    result.push_back(StateVariable::custom(numLevelPops + 2, "Micro turbulence", "velocity"));
    result.push_back(StateVariable::custom(numLevelPops + 3, "temperature", "temperature"));
    for (int p = 1; p <= _numLines; ++p)
        result.push_back(StateVariable::custom(numLevelPops + 3 + p, "radiation field", ""));

    return result;
}

////////////////////////////////////////////////////////////////////

void MolecularLineGasMix::initializeSpecificState(MaterialState* state, double /*metallicity*/, double temperature,
                                                  const Array& params) const
{
    // leave the properties untouched if the cell does not contain any material for this component
    if (state->numberDensity() > 0.)
    {
        // initialization of the convergence of the level populations
        state->setCustom(numLevelPops + 1, 1e99);

        // initialization of microtubulence [m/s]
        state->setCustom(numLevelPops + 2, params.size() ? params[1] : defaultMicroTurbulenceVelocity());

        // if no value was imported, use default value
        state->setCustom(0, params.size() ? params[0] : state->numberDensity() * defaultCollisionPartnerRatios()[0]);

        // if no value was imported, use default value
        // make sure the temperature is at least the local universe CMB temperature
        double effectiveTemperature;

        effectiveTemperature = 20. + (pow(state->custom(numLevelPops + 2), 2.) / Constants::k() * mass());
        state->setTemperature(max(Constants::Tcmb(), effectiveTemperature));
        state->setCustom(numLevelPops + 3, 20.);

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
            if (statePop <= 0.0 && statePop >= state->numberDensity())
                throw FATALERROR("Failed initialization of level population: " + StringUtils::toString(statePop));
            state->setCustom(p + 1, statePop);
        }
    }
}

////////////////////////////////////////////////////////////////////

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
        if (isnan(Cul))
            throw FATALERROR("Failed calculation for a collisional de-excitatipon coefficient Cul "
                             + StringUtils::toString(indexTem));
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
        if (isnan(Clu)) throw FATALERROR("Failed calculation for a collisional excitation coefficeint Clu");
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

        // for debug
        for (int i = 0; i < numLevelPops; ++i)
        {
            for (int j = 0; j < numLevelPops + 1; ++j)
                if (isnan(matrixPop[i][j])) throw FATALERROR("nan after forwading LU decomposition");
        }

        // backward substitution
        for (int i = static_cast<int>(rowSize) - 1; i >= 0; --i)
        {
            for (size_t j = i + 1; j < collumnSize - 1; ++j) solution[i] -= matrixPop[i][j] * solution[j];
            solution[i] /= matrixPop[i][i];
        }

        // for debug
        for (size_t i = 0; i < rowSize; ++i)
            if (isnan(solution[i])) throw FATALERROR("nan after backwading LU decomposition");

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

////////////////////////////////////////////////////////////////////

UpdateStatus MolecularLineGasMix::updateSpecificState(MaterialState* state, const Array& Jv) const
{
    UpdateStatus status;

    // leave the properties untouched if the cell does not contain any material for this component
    if (state->numberDensity() > 0.0)
    {
        auto config = find<Configuration>();

        double numHydroMol = state->custom(0);
        if (isnan(numHydroMol)) throw FATALERROR("Hydrogen is nan");

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
            if (isnan(Jv[_indexLines[i]]))
                throw FATALERROR("Jv is " + StringUtils::toString(Jv[_indexLines[i] - 1]) + " "
                                 + StringUtils::toString(_indexLines[i]) + "/"
                                 + StringUtils::toString(static_cast<int>(Jv.size())));
            else
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

               /* if (integralsum > 1.01 or integralsum < 0.99)
                    find<Log>()->warning(
                        "Integrated Gaussian = " + StringUtils::toString(integralsum)
                        + ", please adjust lineRange in updateSpecificState, wavelengthgrids, or temperature");*/
            }
        }

        for (int i = 0; i < _numLines; ++i)
            if (rf[i] < 0.0)
                throw FATALERROR("Jv of line " + StringUtils::toString(i)
                                 + "is negative: " + StringUtils::toString(Jv[i]));

        // statistical equilibrium matrix for CO level population
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

        // debug: confirm whether there is nan
        for (int i = 0; i < numLevelPops; ++i)
        {
            for (int j = 0; j < numLevelPops + 1; ++j)
                if (isnan(matrixPop[i][j])) throw FATALERROR("nan before gauss elimination 1");
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

        // add the normalization of CO number density in the last row of the matrix.
        for (int j = 0; j <= numLevelPops; ++j)
        {
            if (j == numLevelPops)
                matrixPop[numLevelPops - 1][j] = state->numberDensity();
            else
                matrixPop[numLevelPops - 1][j] = 1.;
        }

        // debug: confirm whether there is nan
        for (int i = 0; i < numLevelPops; ++i)
        {
            for (int j = 0; j < numLevelPops + 1; ++j)
            {
                if (isnan(matrixPop[i][j])) throw FATALERROR("nan before gauss elimination 2");
                //               log->info(StringUtils::toString(i) + ", " + StringUtils::toString(j) + ": "
                // + StringUtils::toString(matrixPop[i][j]));
            }
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

        // for debug: confirm CO number density
        if (total_num >= state->numberDensity() * (1.0 + 1e-5) || total_num <= state->numberDensity() * (1.0 - 1e-5))
            throw FATALERROR("Failed Update total: " + StringUtils::toString(state->numberDensity())
                             + "; after: " + StringUtils::toString(total_num));
        // confirm the convergence of level population
        double converge_value = 0.;
        for (int p = 1; p <= numLevelPops; ++p)
            if (p < numLevelPops) converge_value += abs((levelpop_before[p - 1] / state->custom(p) - 1));
        converge_value /= (numLevelPops - 1);
        state->setCustom(numLevelPops + 1, converge_value);

        // memorize the radiation fields at line centers
        for (int p = 1; p <= _numLines; ++p) state->setCustom(numLevelPops + 3 + p, rf[p - 1]);

        if (converge_value > maxChangeInLevelPopulations()) status.updateNotConverged();
        else status.updateConverged();
    }

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
    return Constants::Mproton() * molecularWeight;
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
        int line_num = -1;

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
                    if (transrate < 0.0) line_num = p;
                }
            }
        }

        // if the opacity is negative, the level population can be inverse
        // (this is not an error but we could not treat the negative opacity)
        if (kappa < 0.0)
            find<Log>()->info("Cell Index " + StringUtils::toString(state->cellIndex())
                              + " kappa= " + StringUtils::toString(kappa) + " Line " + StringUtils::toString(line_num));
    }
    if (kappa < 0.0) kappa = 0.0;

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
    Array centers(_numLines);

    for (int p = 0; p < _numLines; ++p)
        centers[p] = Constants::c()
                     / ((_energyStates[_indexUpRad[_indexRadTrans[p]]] - _energyStates[_indexLowRad[_indexRadTrans[p]]])
                        / Constants::h());
    return centers;
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
            if (isnan(luminosities[p])) throw FATALERROR("line luminosity is nan ");
            if (luminosities[p] < 0.) throw FATALERROR("line luminosity is negaive ");
        }
    }
    return luminosities;
}

////////////////////////////////////////////////////////////////////

double MolecularLineGasMix::indicativeTemperature(const MaterialState* state, const Array& /*Jv*/) const
{
    return state->temperature();
}

////////////////////////////////////////////////////////////////////
