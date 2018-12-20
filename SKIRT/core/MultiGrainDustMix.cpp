/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "MultiGrainDustMix.hpp"
#include "Configuration.hpp"
#include "Constants.hpp"
#include "DisjointWavelengthGrid.hpp"
#include "FatalError.hpp"
#include "GrainComposition.hpp"
#include "GrainSizeDistribution.hpp"
#include "Log.hpp"
#include "NR.hpp"
#include "Parallel.hpp"
#include "ParallelFactory.hpp"
#include "PlanckFunction.hpp"
#include "ProcessManager.hpp"
#include "StoredTable.hpp"
#include "WavelengthGrid.hpp"

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

double MultiGrainDustMix::getOpticalProperties(const Array& lambdav, const Array& thetav,
                                               Array& sigmaabsv, Array& sigmascav, Array& asymmparv,
                                               Table<2>& S11vv, Table<2>& S12vv, Table<2>& S33vv, Table<2>& S34vv)
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
    if (mode == ScatteringMode::HenyeyGreenstein && numMueller!=0)
        throw FATALERROR("HenyeyGreenstein scattering mode prohibits grain populations to offer a Mueller matrix");
    if (mode == ScatteringMode::MaterialPhaseFunction && numMueller!=numPops)
        throw FATALERROR("MaterialPhaseFunction scattering mode requires grain populations to offer a Mueller matrix");
    if (mode == ScatteringMode::SphericalPolarization && numMueller!=numPops)
        throw FATALERROR("SphericalPolarization scattering mode requires grain populations to offer a Mueller matrix");

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
        int numSizes = max(3., 100 * log10(amax/amin));
        Array av(numSizes);      // "a" for each point
        Array dav(numSizes);     // "da" for each point
        Array dndav(numSizes);   // "dnda" for each point
        Array weightv(numSizes); // integration weight for each point (1/2 or 1 in addition to normalization)
        {
            double logamin = log10(amin);
            double logamax = log10(amax);
            double dloga = (logamax-logamin)/(numSizes-1);
            for (int i=0; i!=numSizes; ++i)
            {
                av[i] = pow(10, logamin + i*dloga);
                dav[i] = av[i] * M_LN10 * dloga;
                dndav[i] = population->sizeDistribution()->dnda(av[i]);
                weightv[i] = 1.;
            }
            weightv[0] = weightv[numSizes-1] = 0.5;
        }

        // calculate the mass per hydrogen atom for this population according to the bare size distribution
        // (i.e. without applying any normalization)
        double baremupop = 0.;
        for (int i=0; i!=numSizes; ++i)
        {
            double volume = 4.0*M_PI/3.0 * av[i] * av[i] * av[i];
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
        if (!mupop) throw FATALERROR("Dust grain population of type " + population->composition()->name()
                                  + " has zero dust mass");
        _mupopv.push_back(mupop);
        mu += mupop;

        // remember the size distribution normalization factor for this population
        // and adjust the integration weight for further calculations
        _normv.push_back(mupop/baremupop);
        weightv *= mupop/baremupop;

        // open the stored tables for the basic optical properties
        string opticalPropsName = population->composition()->resourceNameForOpticalProps();
        StoredTable<2> Qabs(this, opticalPropsName, "a(m),lambda(m)", "Qabs(1)");
        StoredTable<2> Qsca(this, opticalPropsName, "a(m),lambda(m)", "Qsca(1)");
        StoredTable<2> g(this, opticalPropsName, "a(m),lambda(m)", "g(1)");

        // calculate the optical properties for each wavelength, and add them to the global total
        for (int ell=0; ell!=numLambda; ++ell)
        {
            double lamdba = lambdav[ell];

            double sumsigmaabs = 0.;
            double sumsigmasca = 0.;
            double sumgsigmasca = 0.;
            for (int i=0; i!=numSizes; ++i)
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

        // handle Mueller matrix coefficients if needed
        if (mode == ScatteringMode::MaterialPhaseFunction || mode == ScatteringMode::SphericalPolarization)
        {
            // open the stored tables for the Mueller matrix coefficients
            string muellerName = population->composition()->resourceNameForMuellerMatrix();
            StoredTable<3> S11, S12, S33, S34;
            S11.open(this, muellerName, "a(m),lambda(m),theta(rad)", "S11(1)");
            if (mode == ScatteringMode::SphericalPolarization)
            {
                S12.open(this, muellerName, "a(m),lambda(m),theta(rad)", "S12(1)");
                S33.open(this, muellerName, "a(m),lambda(m),theta(rad)", "S33(1)");
                S34.open(this, muellerName, "a(m),lambda(m),theta(rad)", "S34(1)");
            }

            // calculate the Mueller matrix coefficients on the requested wavelength and scattering angle grid;
            // this is time-consuming, so we do this in parallel,
            // which means we need to synchronize the resulting tables after accumulation for all populations
            auto log = find<Log>();
            log->info("Integrating Mueller matrix coefficients over the grain size distribution...");
            log->infoSetElapsed(numLambda);
            find<ParallelFactory>()->parallelDistributed()->call(numLambda,
                [&lambdav,&thetav,&av,&dav,&dndav,&weightv,&S11,&S12,&S33,&S34,&S11vv,&S12vv,&S33vv,&S34vv,mode,log]
                (size_t firstIndex, size_t numIndices)
            {
                int numTheta = thetav.size();
                int numSizes = av.size();

                const size_t logProgressChunkSize = 50;
                while (numIndices)
                {
                    size_t currentChunkSize = min(logProgressChunkSize, numIndices);
                    for (size_t ell=firstIndex; ell!=firstIndex+currentChunkSize; ++ell)
                    {
                        double lambda = lambdav[ell];
                        for (int t=0; t!=numTheta; ++t)
                        {
                            double theta = thetav[t];

                            // integrate over the size distribution
                            double sumS11 = 0.;
                            double sumS12 = 0.;
                            double sumS33 = 0.;
                            double sumS34 = 0.;
                            for (int i=0; i!=numSizes; ++i)
                            {
                                double factor = weightv[i] * dndav[i] * dav[i];
                                sumS11 += factor * S11(av[i],lambda,theta);
                                if (mode == ScatteringMode::SphericalPolarization)
                                {
                                    sumS12 += factor * S12(av[i],lambda,theta);
                                    sumS33 += factor * S33(av[i],lambda,theta);
                                    sumS34 += factor * S34(av[i],lambda,theta);
                                }
                            }

                            // accumulate for all populations
                            S11vv(ell,t) += sumS11;
                            if (mode == ScatteringMode::SphericalPolarization)
                            {
                                S12vv(ell,t) += sumS12;
                                S33vv(ell,t) += sumS33;
                                S34vv(ell,t) += sumS34;
                            }
                        }
                    }
                    log->infoIfElapsed("Calculated Mueller matrix coefficients: ", currentChunkSize);
                    firstIndex += currentChunkSize;
                    numIndices -= currentChunkSize;
                }
            });
        }
    }

    // calculate g = gsigmasca / sigmasca
    for (int ell=0; ell!=numLambda; ++ell)
    {
        asymmparv[ell] = sigmascav[ell] ? asymmparv[ell]/sigmascav[ell] : 0.;
    }

    // synchronize the Mueller matrix coefficients between processes, if applicable
    ProcessManager::sumToAll(S11vv.data());
    ProcessManager::sumToAll(S12vv.data());
    ProcessManager::sumToAll(S33vv.data());
    ProcessManager::sumToAll(S34vv.data());

    return mu;
}

////////////////////////////////////////////////////////////////////

size_t MultiGrainDustMix::initializeExtraProperties(const Array& lambdav)
{
    // determine which type(s) of emission we need to support
    auto config = find<Configuration>();
    bool multigrain = config->hasDustEmission();
    bool stochastic = config->hasStochasticDustEmission();

    // perform only if extra properties are required
    if (multigrain)
    {
        // get the number of wavelength grid points
        int numLambda = lambdav.size();

        // determine the total number of size bins
        int numBins = 0;
        for (auto population : _populations) numBins += population->numSizes();

        // resize the arrays needed for stochastic emissivity
        // (the equilibrium temperature calculator does not need to know the number of bins in advance)
        _sigmaabsvv.resize(numBins, numLambda);
        if (stochastic)
        {
            _btocv.resize(numBins);
            _massv.resize(numBins);
            _enthalpyv.resize(_populations.size());
        }

        // loop over all populations and process size bins for each
        int c = 0;      // population index
        int b = 0;      // running size bin index
        for (auto population : _populations)
        {
            // open the stored table for the absorption cross sections
            string opticalPropsName = population->composition()->resourceNameForOpticalProps();
            StoredTable<2> Qabs(this, opticalPropsName, "a(m),lambda(m)", "Qabs(1)");

            // construct the size bins (i.e. the bin border points) for this population
            int numPopBins = population->numSizes();
            double amin = population->sizeDistribution()->amin();
            double amax = population->sizeDistribution()->amax();
            Array aborderv;
            NR::buildLogGrid(aborderv, amin, amax, numPopBins);

            // loop over the size bins for this population
            for (int bb = 0; bb!=numPopBins; ++bb)
            {
                // create an integration grid over grain size within this bin
                int numSizes = max(3., 100 * log10(amax/amin));
                Array av(numSizes);      // "a" for each point
                Array dav(numSizes);     // "da" for each point
                Array dndav(numSizes);   // "dnda" for each point
                Array weightv(numSizes); // integration weight for each point (1/2 or 1 in addition to normalization)
                {
                    double logamin = log10(aborderv[bb]);
                    double logamax = log10(aborderv[bb+1]);
                    double dloga = (logamax-logamin)/(numSizes-1);
                    for (int i=0; i!=numSizes; ++i)
                    {
                        av[i] = pow(10, logamin + i*dloga);
                        dav[i] = av[i] * M_LN10 * dloga;
                        dndav[i] = population->sizeDistribution()->dnda(av[i]);
                        weightv[i] = _normv[c];
                    }
                    weightv[0] *= 0.5;
                    weightv[numSizes-1] *= 0.5;
                }

                // size-integrate the absorption cross sections for this bin
                for (int ell=0; ell!=numLambda; ++ell)
                {
                    double lamdba = lambdav[ell];
                    double sum = 0.;
                    for (int i=0; i!=numSizes; ++i)
                    {
                        double area = M_PI * av[i] * av[i];
                        sum += weightv[i] * dndav[i] * area * Qabs(av[i], lamdba) * dav[i];
                    }
                    _sigmaabsvv(b,ell) = sum;
                }

                // setup the equilibrium temperature calculator for this bin
                _tempCalcv.precalculate(this, lambdav, _sigmaabsvv[b]);

                // if needed, also precalculate the info for stochastic emissivity calculations
                if (stochastic)
                {
                    // remember the mapping from bins to populations
                    _btocv[b] = c;

                    // calculate the mean grain mass for this bin
                    double sum1 = 0.;
                    double sum2 = 0.;
                    for (int i=0; i!=numSizes; ++i)
                    {
                        double volume = 4.0*M_PI/3.0 * av[i] * av[i] * av[i];
                        sum1 += weightv[i] * dndav[i] * volume * dav[i];
                        sum2 += weightv[i] * dndav[i] * dav[i];
                    }
                    _massv[b] = population->composition()->bulkDensity() * sum1/sum2;
                }

                // increment the running bin index
                b++;
            }

            if (stochastic)
            {
                // open the enthalpy stored table for this population
                string enthalpyName = population->composition()->resourceNameForEnthalpies();
                _enthalpyv[c] = new StoredTable<1>(this, enthalpyName, "T(K)", "h(J/m3)");
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
    allocatedBytes += _tempCalcv.allocatedBytes();
    allocatedBytes += _sigmaabsvv.size() * sizeof(double);
    allocatedBytes += _btocv.size() * sizeof(_btocv[0]);
    allocatedBytes += _massv.size() * sizeof(_massv[0]);
    allocatedBytes += _enthalpyv.size() * sizeof(_enthalpyv[0]);
    return allocatedBytes;
}

////////////////////////////////////////////////////////////////////

MultiGrainDustMix::~MultiGrainDustMix()
{
    for (auto enthalpy : _enthalpyv) delete enthalpy;
}

////////////////////////////////////////////////////////////////////

Array MultiGrainDustMix::emissivity(const Array& Jv) const
{
    // get the output wavelength grid
    auto wavelengthGrid = find<Configuration>()->dustEmissionWLG();
    int numWavelengths = wavelengthGrid->numBins();

    // calculate the black-body emissivity spectrum
    Array ev(numWavelengths);
    int numBins = _tempCalcv.numBins();
    for (int b=0; b!=numBins; ++b)
    {
        double T = _tempCalcv.equilibriumTemperature(b, Jv);
        PlanckFunction B(T);
        for (int ell=0; ell<numWavelengths; ell++)
        {
            double lambda = wavelengthGrid->wavelength(ell);
            ev[ell] += _sigmaabsvv(b, indexForLambda(lambda)) * B(lambda);
        }
    }
    return ev;
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

double MultiGrainDustMix::populationMass(int c) const
{
    return _mupopv[c];
}

/*
////////////////////////////////////////////////////////////////////

int MultiGrainDustMix::numBins() const
{
    return _tempCalcv.numBins();
}

////////////////////////////////////////////////////////////////////

double MultiGrainDustMix::binEquilibriumTemperature(int b, const Array& Jv) const
{
    return _tempCalcv.equilibriumTemperature(b, Jv);
}

////////////////////////////////////////////////////////////////////

double MultiGrainDustMix::binSectionAbs(int b, double lambda) const
{
    return _sigmaabsvv(b, indexForLambda(lambda));
}

////////////////////////////////////////////////////////////////////

string MultiGrainDustMix::binGrainType(int b) const
{
    int c = _btocv[b];
    return _populations[c]->composition()->name();
}

////////////////////////////////////////////////////////////////////

double MultiGrainDustMix::binMeanMass(int b) const
{
    return _massv[b];
}

////////////////////////////////////////////////////////////////////

double MultiGrainDustMix::binEnthalpy(int b, double T) const
{
    int c = _btocv[b];
    return _massv[b] * _enthalpyv[c]->operator()(T) / _populations[c]->composition()->bulkDensity();
}

////////////////////////////////////////////////////////////////////

double MultiGrainDustMix::maxEnthalpyTemperature() const
{
    Range range(0., 3.);
    for (auto enthalpy : _enthalpyv) range.extend(enthalpy->axisRange<0>());
    return range.max();
}
*/
////////////////////////////////////////////////////////////////////
