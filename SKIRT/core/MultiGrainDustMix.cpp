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

double MultiGrainDustMix::getOpticalProperties(const Array& lambdav, const Array& thetav,
                                               Array& sigmaabsv, Array& sigmascav, Array& asymmparv,
                                               Table<2>& S11vv, Table<2>& S12vv, Table<2>& S33vv, Table<2>& S34vv) const
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
    int numTheta = thetav.size();

    // dust mass per hydrogen atom accumulated over all populations
    double mu = 0.;

    // accumulate the relevant properties over all populations
    for (auto population : _populations)
    {
        // construct a grain size integration grid for this population
        double amin = population->sizeDistribution()->amin();
        double amax = population->sizeDistribution()->amax();
        int numSizes = max(3., 200 * log10(amax/amin));
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
        mu += mupop;

        // adjust the integration weight for further calculations by the normalization factor
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

            // calculate the Mueller matrix coefficients on the requested wavelength and scattering angle grid
            for (int ell=0; ell!=numLambda; ++ell)
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
        }
    }

    // calculate g = gsigmasca / sigmasca
    for (int ell=0; ell!=numLambda; ++ell)
    {
        asymmparv[ell] = sigmascav[ell] ? asymmparv[ell]/sigmascav[ell] : 0.;
    }

    return mu;
}

////////////////////////////////////////////////////////////////////
