/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "MultiGrainDustMix.hpp"
#include "GrainPopulation.hpp"
#include "NR.hpp"

////////////////////////////////////////////////////////////////////

void MultiGrainDustMix::addPopulation(const GrainPopulation* /*population*/)
{
}

////////////////////////////////////////////////////////////////////

void MultiGrainDustMix::addPopulation(GrainComposition* composition, GrainSizeDistribution* sizeDistribution,
                                      int numSizes, GrainPopulation::NormalizationType normType, double normValue)
{
    addPopulation(new GrainPopulation(this, composition, sizeDistribution, numSizes, normType, normValue));
}

////////////////////////////////////////////////////////////////////

MaterialMix::ScatteringMode MultiGrainDustMix::scatteringMode() const
{
    return MaterialMix::ScatteringMode::HenyeyGreenstein;
}

////////////////////////////////////////////////////////////////////

double MultiGrainDustMix::mass() const
{
    return 0.;
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
