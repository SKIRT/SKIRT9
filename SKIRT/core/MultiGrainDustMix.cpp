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
#include "NR.hpp"
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

void MultiGrainDustMix::setupSelfAfter()
{
    // get the number of grain populations
    int numPops = _populations.size();
    if (!numPops) throw FATALERROR("Dust mix must have at least one grain population");

    // construct a wavelength grid with an appropriate range and 300 wavelengths per decade;
    // this might be optimized by varying the wavelength resolution depending on the subrange
    auto range = find<Configuration>()->simulationWavelengthRange();
    int numWavelengths = max(3., 300 * log10(range.max()/range.min()));
    NR::buildLogGrid(_lambdav, range.min(), range.max(), numWavelengths);

    // resize the arrays for properties tabulated on wavelength
    _sigmaabsv.resize(numWavelengths);
    _sigmascav.resize(numWavelengths);
    _asymmparv.resize(numWavelengths);

    // accumulate the relevant properties over all populations
    for (auto population : _populations)
    {
        // open the stored tables
        string opticalPropsName = population->composition()->resourceNameForOpticalProps();
        StoredTable<2> Qabs(this, opticalPropsName, "lambda(m),a(m)", "Qabs(1)");
        StoredTable<2> Qsca(this, opticalPropsName, "lambda(m),a(m)", "Qsca(1)");
        StoredTable<2> g(this, opticalPropsName, "lambda(m),a(m)", "g(1)");

        // construct a grain size integration grid for this population
        double amin = population->sizeDistribution()->amin();
        double amax = population->sizeDistribution()->amax();
        int numSizes = max(3., 200 * log10(amax/amin));
        Array av(numSizes);      // "a" for each point
        Array dav(numSizes);     // "da" for each point
        Array dndav(numSizes);   // "dnda" for each point
        Array weightv(numSizes); // integration weight for each point (1/2 or 1)
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
        double baremu = 0.;
        for (int i=0; i!=numSizes; ++i)
        {
            double volume = 4.0*M_PI/3.0 * av[i] * av[i] * av[i];
            baremu += weightv[i] * dndav[i] * volume * dav[i];
        }
        baremu *= population->composition()->bulkDensity();

        // determine the actual mass per hydrogen atom for this population after applying normalization,
        // and add it to the global total
        double mu = 0.;
        switch (population->normalizationType())
        {
        case GrainPopulation::NormalizationType::DustMassPerHydrogenAtom:
            mu = population->dustMassPerHydrogenAtom();
            break;
        case GrainPopulation::NormalizationType::DustMassPerHydrogenMass:
            mu = population->dustMassPerHydrogenMass() * Constants::Mproton();
            break;
        case GrainPopulation::NormalizationType::FactorOnSizeDistribution:
            mu = baremu * population->factorOnSizeDistribution();
            break;
        }
        if (!mu) throw FATALERROR("Dust grain population of type " + population->composition()->name()
                                  + " has zero dust mass");
        _mu += mu;

        // adjust the integration weight for further calculations by the normalization factor
        weightv *= mu/baremu;

        // calculate the optical properties for each wavelength, and add them to the global total
        for (int ell=0; ell!=numWavelengths; ++ell)
        {
            double lamdba = _lambdav[ell];

            double sumsigmaabs = 0.;
            double sumsigmasca = 0.;
            double sumgsigmasca = 0.;
            for (int i=0; i!=numSizes; ++i)
            {
                double area = M_PI * av[i] * av[i];
                double factor = weightv[i] * dndav[i] * area * dav[i];
                double sigmaabs = factor * Qabs(lamdba ,av[i]);
                double sigmasca = factor * Qsca(lamdba, av[i]);
                double gsigmasca = sigmasca * g(lamdba, av[i]);
                sumsigmaabs += sigmaabs;
                sumsigmasca += sigmasca;
                sumgsigmasca += gsigmasca;
            }
            _sigmaabsv[ell] += sumsigmaabs;
            _sigmascav[ell] += sumsigmasca;
            _asymmparv[ell] += sumgsigmasca;  // need to divide by sigmasca after accumulation for all populations
        }
    }

    // calculate g = gsigmasca / sigmasca
    for (int ell=0; ell!=numWavelengths; ++ell)
    {
        _asymmparv[ell] = _sigmascav[ell] ? _asymmparv[ell]/_sigmascav[ell] : 0.;
    }

    // TO DO
    DustMix::setupSelfAfter();
}

////////////////////////////////////////////////////////////////////

MaterialMix::ScatteringMode MultiGrainDustMix::scatteringMode() const
{
    return MaterialMix::ScatteringMode::HenyeyGreenstein;
}

////////////////////////////////////////////////////////////////////

double MultiGrainDustMix::mass() const
{
    return _mu;
}

////////////////////////////////////////////////////////////////////

double MultiGrainDustMix::sectionAbsSelf(double lambda) const
{
    return NR::clampedValue<NR::interpolateLogLog>(lambda, _lambdav, _sigmaabsv);
}

////////////////////////////////////////////////////////////////////

double MultiGrainDustMix::sectionScaSelf(double lambda) const
{
    return NR::clampedValue<NR::interpolateLogLog>(lambda, _lambdav, _sigmascav);
}

////////////////////////////////////////////////////////////////////

double MultiGrainDustMix::asymmpar(double lambda) const
{
    return NR::clampedValue<NR::interpolateLogLin>(lambda, _lambdav, _asymmparv);
}

////////////////////////////////////////////////////////////////////

double MultiGrainDustMix::phaseFunctionValueForCosine(double /*lambda*/, double /*costheta*/) const
{
    return 0.;
}

////////////////////////////////////////////////////////////////////

double MultiGrainDustMix::generateCosineFromPhaseFunction(double /*lambda*/) const
{
    return 0.;
}

////////////////////////////////////////////////////////////////////

double MultiGrainDustMix::phaseFunctionValue(double /*lambda*/, double /*theta*/, double /*phi*/, const StokesVector* /*sv*/) const
{
    return 0.;
}

////////////////////////////////////////////////////////////////////

std::pair<double,double> MultiGrainDustMix::generateAnglesFromPhaseFunction(double /*lambda*/, const StokesVector* /*sv*/) const
{
    return std::make_pair(0.,0.);
}

////////////////////////////////////////////////////////////////////

void MultiGrainDustMix::applyMueller(double /*lambda*/, double /*theta*/, StokesVector* /*sv*/) const
{
}

////////////////////////////////////////////////////////////////////
