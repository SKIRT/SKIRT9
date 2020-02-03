/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "MultiGrainDustMix.hpp"
#include "Configuration.hpp"
#include "Constants.hpp"
#include "FatalError.hpp"
#include "GrainComposition.hpp"
#include "GrainSizeDistribution.hpp"
#include "Log.hpp"
#include "NR.hpp"
#include "Parallel.hpp"
#include "ParallelFactory.hpp"
#include "ProcessManager.hpp"
#include "StoredTable.hpp"

////////////////////////////////////////////////////////////////////

void MultiGrainDustMix::addPopulation(const GrainPopulation* population)
{
    _populations.push_back(population);
}

////////////////////////////////////////////////////////////////////

void MultiGrainDustMix::addPopulation(GrainComposition* composition, GrainSizeDistribution* sizeDistribution,
                                      int numSizes, GrainPopulation::NormalizationType normType, double normValue)
{
    addPopulation(new GrainPopulation(this, composition, sizeDistribution, numSizes, normType, normValue));
}

////////////////////////////////////////////////////////////////////

double MultiGrainDustMix::getOpticalProperties(const Array& lambdav, const Array& thetav, Array& sigmaabsv,
                                               Array& sigmascav, Array& asymmparv, Table<2>& S11vv, Table<2>& S12vv,
                                               Table<2>& S33vv, Table<2>& S34vv, ArrayTable<2>& sigmaabsvv,
                                               ArrayTable<2>& sigmaabspolvv)
{
    // get the scattering mode advertised by this dust mix
    auto mode = scatteringMode();

    // get the number of grain populations
    int numPops = _populations.size();
    if (!numPops) throw FATALERROR("Dust mix must have at least one grain population");

    // count the number of populations that offer a Mueller matrix
    int numMueller = 0;
    for (auto population : _populations)
        if (!population->composition()->resourceNameForMuellerMatrix().empty()) numMueller++;

    // verify the Mueller matrix support by all grain populations
    if (mode == ScatteringMode::HenyeyGreenstein && numMueller != 0)
        throw FATALERROR("HenyeyGreenstein scattering mode prohibits grain populations to offer a Mueller matrix");
    if (mode == ScatteringMode::MaterialPhaseFunction && numMueller != numPops)
        throw FATALERROR("MaterialPhaseFunction scattering mode requires grain populations to offer a Mueller matrix");
    if ((mode == ScatteringMode::SphericalPolarization || mode == ScatteringMode::SpheroidalPolarization)
        && numMueller != numPops)
        throw FATALERROR("SphericalPolarization scattering mode requires grain populations to offer a Mueller matrix");

    // define shortcuts for various mode checks
    bool needSpheroidalPolarization = (mode == ScatteringMode::SpheroidalPolarization);
    bool needSphericalPolarization = needSpheroidalPolarization || (mode == ScatteringMode::SphericalPolarization);
    bool needMaterialPhaseFunction = needSphericalPolarization || (mode == ScatteringMode::MaterialPhaseFunction);

    // get the number of requested grid points
    int numLambda = lambdav.size();

    // dust mass per hydrogen atom accumulated over all populations
    double mu = 0.;

    // accumulate the relevant properties over all populations
    for (auto population : _populations)
    {
        // construct a grain size integration grid for this population
        double amin = population->sizeDistribution()->amin();
        double amax = population->sizeDistribution()->amax();
        int numSizes = max(3., 100 * log10(amax / amin));
        Array av(numSizes);       // "a" for each point
        Array dav(numSizes);      // "da" for each point
        Array dndav(numSizes);    // "dnda" for each point
        Array weightv(numSizes);  // integration weight for each point (1/2 or 1 in addition to normalization)
        {
            double logamin = log10(amin);
            double logamax = log10(amax);
            double dloga = (logamax - logamin) / (numSizes - 1);
            for (int i = 0; i != numSizes; ++i)
            {
                av[i] = pow(10, logamin + i * dloga);
                dav[i] = av[i] * M_LN10 * dloga;
                dndav[i] = population->sizeDistribution()->dnda(av[i]);
                weightv[i] = 1.;
            }
            weightv[0] = weightv[numSizes - 1] = 0.5;
        }

        // calculate the mass per hydrogen atom for this population according to the bare size distribution
        // (i.e. without applying any normalization)
        double baremupop = 0.;
        for (int i = 0; i != numSizes; ++i)
        {
            double volume = 4.0 * M_PI / 3.0 * av[i] * av[i] * av[i];
            baremupop += weightv[i] * dndav[i] * volume * dav[i];
        }
        baremupop *= population->composition()->bulkDensity();

        // determine the actual mass per hydrogen atom for this population after applying normalization,
        // and add it to the global total
        double mupop = 0.;
        switch (population->normalizationType())
        {
            case GrainPopulation::NormalizationType::DustMassPerHydrogenAtom:
                mupop = population->dustMassPerHydrogenAtom();
                break;
            case GrainPopulation::NormalizationType::DustMassPerHydrogenMass:
                mupop = population->dustMassPerHydrogenMass() * Constants::Mproton();
                break;
            case GrainPopulation::NormalizationType::FactorOnSizeDistribution:
                mupop = baremupop * population->factorOnSizeDistribution();
                break;
        }
        if (!mupop)
            throw FATALERROR("Dust grain population of type " + population->composition()->name()
                             + " has zero dust mass");
        _mupopv.push_back(mupop);
        mu += mupop;

        // remember the size distribution normalization factor for this population
        // and adjust the integration weight for further calculations
        _normv.push_back(mupop / baremupop);
        weightv *= mupop / baremupop;

        // open the stored tables for the basic optical properties
        string opticalPropsName = population->composition()->resourceNameForOpticalProps();
        StoredTable<2> Qabs(this, opticalPropsName, "a(m),lambda(m)", "Qabs(1)");
        StoredTable<2> Qsca(this, opticalPropsName, "a(m),lambda(m)", "Qsca(1)");
        StoredTable<2> g(this, opticalPropsName, "a(m),lambda(m)", "g(1)");

        // calculate the optical properties for each wavelength, and add them to the global total
        for (int ell = 0; ell != numLambda; ++ell)
        {
            double lamdba = lambdav[ell];

            double sumsigmaabs = 0.;
            double sumsigmasca = 0.;
            double sumgsigmasca = 0.;
            for (int i = 0; i != numSizes; ++i)
            {
                double area = M_PI * av[i] * av[i];
                double factor = weightv[i] * dndav[i] * area * dav[i];
                double sigmaabs = factor * Qabs(av[i], lamdba);
                double sigmasca = factor * Qsca(av[i], lamdba);
                double gsigmasca = sigmasca * g(av[i], lamdba);
                sumsigmaabs += sigmaabs;
                sumsigmasca += sigmasca;
                sumgsigmasca += gsigmasca;
            }
            sigmaabsv[ell] += sumsigmaabs;
            sigmascav[ell] += sumsigmasca;
            asymmparv[ell] += sumgsigmasca;  // need to divide by sigmasca after accumulation for all populations
        }

        // handle Mueller matrix coefficients and spheroidal grain properties if needed
        if (needMaterialPhaseFunction)
        {
            // open the stored tables for the Mueller matrix coefficients
            StoredTable<3> S11, S12, S33, S34;
            string muellerName = population->composition()->resourceNameForMuellerMatrix();
            S11.open(this, muellerName, "a(m),lambda(m),theta(rad)", "S11(1)");
            if (needSphericalPolarization)
            {
                S12.open(this, muellerName, "a(m),lambda(m),theta(rad)", "S12(1)");
                S33.open(this, muellerName, "a(m),lambda(m),theta(rad)", "S33(1)");
                S34.open(this, muellerName, "a(m),lambda(m),theta(rad)", "S34(1)");
            }

            // open the stored tables for the spheroidal grain efficiencies, if present for this population
            StoredTable<3> sQabs, sQabspol;
            string spheroidalPropsName = population->composition()->resourceNameForSpheroidalEmission();
            bool hasSpheroidalPolarization = needSpheroidalPolarization && !spheroidalPropsName.empty();
            if (hasSpheroidalPolarization)
            {
                sQabs.open(this, spheroidalPropsName, "a(m),lambda(m),theta(rad)", "Qabs(1)", true, false);
                sQabspol.open(this, spheroidalPropsName, "a(m),lambda(m),theta(rad)", "Qabspol(1)", true, false);
            }

            // calculate the required grain properties on the requested wavelength and scattering angle grid;
            // this is time-consuming, so we do this in parallel,
            // which means we need to synchronize the resulting tables after accumulation for all populations
            auto log = find<Log>();
            log->info("Integrating grain properties over the grain size distribution...");
            log->infoSetElapsed(numLambda);
            find<ParallelFactory>()->parallelDistributed()->call(
                numLambda,
                [log, needSphericalPolarization, needSpheroidalPolarization, hasSpheroidalPolarization, &lambdav,
                 &thetav, &av, &dav, &dndav, &weightv, &S11, &S12, &S33, &S34, &S11vv, &S12vv, &S33vv, &S34vv, &Qabs,
                 &sQabs, &sQabspol, &sigmaabsvv, &sigmaabspolvv](size_t firstIndex, size_t numIndices) {
                    int numTheta = thetav.size();
                    int numSizes = av.size();

                    const size_t logProgressChunkSize = 50;
                    while (numIndices)
                    {
                        size_t currentChunkSize = min(logProgressChunkSize, numIndices);
                        for (size_t ell = firstIndex; ell != firstIndex + currentChunkSize; ++ell)
                        {
                            double lambda = lambdav[ell];
                            for (int t = 0; t != numTheta; ++t)
                            {
                                double theta = thetav[t];

                                // integrate over the size distribution
                                double sumS11 = 0.;
                                double sumS12 = 0.;
                                double sumS33 = 0.;
                                double sumS34 = 0.;
                                double sumabs = 0.;
                                double sumabspol = 0.;
                                for (int i = 0; i != numSizes; ++i)
                                {
                                    double a = av[i];
                                    double factor = weightv[i] * dndav[i] * dav[i];
                                    sumS11 += factor * S11(av[i], lambda, theta);
                                    if (needSphericalPolarization)
                                    {
                                        sumS12 += factor * S12(a, lambda, theta);
                                        sumS33 += factor * S33(a, lambda, theta);
                                        sumS34 += factor * S34(a, lambda, theta);
                                    }
                                    if (needSpheroidalPolarization)
                                    {
                                        factor *= M_PI * a * a;
                                        if (hasSpheroidalPolarization)
                                        {
                                            sumabs += factor * sQabs(a, lambda, theta);
                                            sumabspol += factor * sQabspol(a, lambda, theta);
                                        }
                                        else
                                        {
                                            // for spherical grains:
                                            //    Qabs(a,lambda,theta) = Qabs(a,lambda)
                                            //    Qabspol(a,lambda,theta) = 0
                                            sumabs += factor * Qabs(a, lambda);
                                        }
                                    }
                                }

                                // accumulate for all populations
                                S11vv(ell, t) += sumS11;
                                if (needSphericalPolarization)
                                {
                                    S12vv(ell, t) += sumS12;
                                    S33vv(ell, t) += sumS33;
                                    S34vv(ell, t) += sumS34;
                                }
                                if (needSpheroidalPolarization)
                                {
                                    sigmaabsvv(ell, t) += sumabs;
                                    sigmaabspolvv(ell, t) += sumabspol;
                                }
                            }
                        }
                        log->infoIfElapsed("Calculated grain properties: ", currentChunkSize);
                        firstIndex += currentChunkSize;
                        numIndices -= currentChunkSize;
                    }
                });
        }
    }

    // calculate g = gsigmasca / sigmasca
    for (int ell = 0; ell != numLambda; ++ell)
    {
        asymmparv[ell] = sigmascav[ell] ? asymmparv[ell] / sigmascav[ell] : 0.;
    }

    // synchronize the Mueller matrix coefficients between processes, if applicable
    if (needMaterialPhaseFunction)
    {
        ProcessManager::sumToAll(S11vv.data());
        if (needSphericalPolarization)
        {
            ProcessManager::sumToAll(S12vv.data());
            ProcessManager::sumToAll(S33vv.data());
            ProcessManager::sumToAll(S34vv.data());
        }
    }

    // synchronize the spheroidal grain cross sections between processes, if applicable
    if (needSpheroidalPolarization)
    {
        for (int ell = 0; ell != numLambda; ++ell)
        {
            ProcessManager::sumToAll(sigmaabsvv[ell]);
            ProcessManager::sumToAll(sigmaabspolvv[ell]);
        }
    }
    return mu;
}

////////////////////////////////////////////////////////////////////

size_t MultiGrainDustMix::initializeExtraProperties(const Array& lambdav)
{
    // determine which type(s) of emission we need to support
    auto config = find<Configuration>();
    _multigrain = config->hasDustEmission();
    _stochastic = config->hasStochasticDustEmission();

    // perform only if extra properties are required
    if (_multigrain)
    {
        // get the number of wavelength grid points
        size_t numLambda = lambdav.size();

        // allocate temporary array for size-bin-integrated absorption cross sections
        Array sigmaabsv(numLambda);

        // loop over all populations and process size bins for each
        int c = 0;  // population index
        int b = 0;  // running size bin index
        for (auto population : _populations)
        {
            // open the absorption cross section stored table for this population
            string opticalPropsName = population->composition()->resourceNameForOpticalProps();
            StoredTable<2> Qabs(this, opticalPropsName, "a(m),lambda(m)", "Qabs(1)");

            // if applicable, open the enthalpy stored table for this population
            StoredTable<1> enthalpy;
            if (_stochastic)
            {
                string enthalpyName = population->composition()->resourceNameForEnthalpies();
                enthalpy.open(this, enthalpyName, "T(K)", "h(J/m3)");
            }

            // construct the size bins (i.e. the bin border points) for this population
            int numPopBins = population->numSizes();
            double amin = population->sizeDistribution()->amin();
            double amax = population->sizeDistribution()->amax();
            Array aborderv;
            NR::buildLogGrid(aborderv, amin, amax, numPopBins);

            // loop over the size bins for this population
            for (int bb = 0; bb != numPopBins; ++bb)
            {
                // create an integration grid over grain size within this bin
                int numSizes = max(3., 100 * log10(amax / amin));
                Array av(numSizes);       // "a" for each point
                Array dav(numSizes);      // "da" for each point
                Array dndav(numSizes);    // "dnda" for each point
                Array weightv(numSizes);  // integration weight for each point (1/2 or 1 in addition to normalization)
                {
                    double logamin = log10(aborderv[bb]);
                    double logamax = log10(aborderv[bb + 1]);
                    double dloga = (logamax - logamin) / (numSizes - 1);
                    for (int i = 0; i != numSizes; ++i)
                    {
                        av[i] = pow(10, logamin + i * dloga);
                        dav[i] = av[i] * M_LN10 * dloga;
                        dndav[i] = population->sizeDistribution()->dnda(av[i]);
                        weightv[i] = _normv[c];
                    }
                    weightv[0] *= 0.5;
                    weightv[numSizes - 1] *= 0.5;
                }

                // size-integrate the absorption cross sections for this bin
                // this can take a few seconds for all populations/size bins combined,
                // so we parallelize the loop but there is no reason to log progress
                sigmaabsv = 0;  // clear array in case calculation is distributed over multiple processes
                find<ParallelFactory>()->parallelDistributed()->call(
                    numLambda,
                    [&lambdav, &av, &dav, &dndav, &weightv, &Qabs, &sigmaabsv](size_t firstIndex, size_t numIndices) {
                        size_t numSizes = av.size();
                        for (size_t ell = firstIndex; ell != firstIndex + numIndices; ++ell)
                        {
                            double sum = 0.;
                            for (size_t i = 0; i != numSizes; ++i)
                            {
                                double area = M_PI * av[i] * av[i];
                                sum += weightv[i] * dndav[i] * area * Qabs(av[i], lambdav[ell]) * dav[i];
                            }
                            sigmaabsv[ell] = sum;
                        }
                    });
                ProcessManager::sumToAll(sigmaabsv);

                // setup the appropriate emissivity calculator for this bin
                if (_stochastic)
                {
                    // calculate the mean grain mass for this bin
                    double sum1 = 0.;
                    double sum2 = 0.;
                    for (int i = 0; i != numSizes; ++i)
                    {
                        double volume = 4.0 * M_PI / 3.0 * av[i] * av[i] * av[i];
                        sum1 += weightv[i] * dndav[i] * volume * dav[i];
                        sum2 += weightv[i] * dndav[i] * dav[i];
                    }
                    double bulkDensity = population->composition()->bulkDensity();
                    double meanMass = sum2 ? bulkDensity * sum1 / sum2 : 0.;

                    // get the grain type for this population
                    string grainType = population->composition()->name();

                    // setup the calculator for this bin
                    _calcSt.precalculate(this, lambdav, sigmaabsv, grainType, bulkDensity, meanMass, enthalpy);
                }
                else
                {
                    _calcEq.precalculate(this, lambdav, sigmaabsv);
                }

                // increment the running bin index
                b++;
            }

            // increment the population index
            c++;
        }
    }

    // determine the allocated number of bytes
    size_t allocatedBytes = 0;
    allocatedBytes += _populations.size() * sizeof(_populations[0]);
    allocatedBytes += _mupopv.size() * sizeof(_mupopv[0]);
    allocatedBytes += _normv.size() * sizeof(_normv[0]);
    allocatedBytes += _calcEq.allocatedBytes();
    allocatedBytes += _calcSt.allocatedBytes();
    return allocatedBytes;
}

////////////////////////////////////////////////////////////////////

Array MultiGrainDustMix::emissivity(const Array& Jv) const
{
    // use the appropriate emissivity calculator
    if (_stochastic)
        return _calcSt.emissivity(Jv);
    else
        return _calcEq.emissivity(Jv);
}

////////////////////////////////////////////////////////////////////

int MultiGrainDustMix::numPopulations() const
{
    return _populations.size();
}

////////////////////////////////////////////////////////////////////

string MultiGrainDustMix::populationGrainType(int c) const
{
    return _populations[c]->composition()->name();
}

////////////////////////////////////////////////////////////////////

Range MultiGrainDustMix::populationSizeRange(int c) const
{
    auto sd = _populations[c]->sizeDistribution();
    return Range(sd->amin(), sd->amax());
}

////////////////////////////////////////////////////////////////////

const GrainSizeDistribution* MultiGrainDustMix::populationSizeDistribution(int c) const
{
    return _populations[c]->sizeDistribution();
}

////////////////////////////////////////////////////////////////////

double MultiGrainDustMix::populationMass(int c) const
{
    return _mupopv[c];
}

////////////////////////////////////////////////////////////////////

void MultiGrainDustMix::getSizeBinProperties(const Array& lambdav, vector<Array>& sizevv,
                                             vector<Array>& numberDensityFractionvv,
                                             vector<vector<Array>>& qabsvvv) const
{
    // get the number of grain populations
    int numPops = _populations.size();
    if (!numPops) throw FATALERROR("Dust mix must have at least one grain population");

    // get the number of requested grid points
    int numLambda = lambdav.size();

    // total number density over all populations
    double mixnumberdens = 0;

    // loop over all populations and process size bins for each (this was copy-pasted from
    // initializeExtraProperties(), and then edited down)
    int c = 0;  // population index
    for (auto population : _populations)
    {
        // open the absorption cross section stored table for this population
        string opticalPropsName = population->composition()->resourceNameForOpticalProps();
        StoredTable<2> Qabs(this, opticalPropsName, "a(m),lambda(m)", "Qabs(1)");

        // construct the size bins (i.e. the bin border points) for this population
        int numPopBins = population->numSizes();
        double amin = population->sizeDistribution()->amin();
        double amax = population->sizeDistribution()->amax();
        Array aborderv;
        NR::buildLogGrid(aborderv, amin, amax, numPopBins);

        // Allocate space for the results
        sizevv.push_back(Array(numPopBins));                                  // indexed on [c][b]
        numberDensityFractionvv.push_back(Array(numPopBins));                 // indexed on [c][b]
        qabsvvv.push_back(std::vector<Array>(numPopBins, Array(numLambda)));  // indexed on [c][b][ell]

        // loop over the size bins for this population
        for (int bb = 0; bb != numPopBins; ++bb)
        {
            // create an integration grid over grain size within this bin
            int numSizes = max(3., 100 * log10(amax / amin));
            Array av(numSizes);       // "a" for each point
            Array dav(numSizes);      // "da" for each point
            Array dndav(numSizes);    // "dnda" for each point
            Array weightv(numSizes);  // integration weight for each point (1/2 or 1 in addition to normalization)
            {
                double logamin = log10(aborderv[bb]);
                double logamax = log10(aborderv[bb + 1]);
                double dloga = (logamax - logamin) / (numSizes - 1);
                for (int i = 0; i != numSizes; ++i)
                {
                    av[i] = pow(10, logamin + i * dloga);
                    dav[i] = av[i] * M_LN10 * dloga;
                    dndav[i] = population->sizeDistribution()->dnda(av[i]);
                    weightv[i] = _normv[c];
                }
                weightv[0] *= 0.5;
                weightv[numSizes - 1] *= 0.5;
            }

            // Integrate the average size and the total number density for this bin.
            double aSum = 0.;
            double n = 0.;
            for (int i = 0; i != numSizes; ++i)
            {
                aSum += weightv[i] * dndav[i] * av[i] * dav[i];
                n += weightv[i] * dndav[i] * dav[i];
            }
            sizevv[c][bb] = aSum / n;
            numberDensityFractionvv[c][bb] = n;
            mixnumberdens += n;

            // Integrate Qabs for this bin. TODO: parallelize as shown in
            // initializeExtraProperties if slow.
            for (int ell = 0; ell != numLambda; ++ell)
            {
                double sum = 0.;
                for (int i = 0; i != numSizes; ++i)
                {
                    sum += weightv[i] * dndav[i] * Qabs(av[i], lambdav[ell]) * dav[i];
                }
                qabsvvv[c][bb][ell] = sum;
            }
        }

        // increment the population index
        c++;
    }

    // convert the stored number densities to number density fractions relative to the entire mix
    for (Array& numberDensityFractionv : numberDensityFractionvv) numberDensityFractionv /= mixnumberdens;
}

////////////////////////////////////////////////////////////////////
