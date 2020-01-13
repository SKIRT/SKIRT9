/*//////////////////////////////////////////////////////////////////
////       SKIRT -- an advanced radiative transfer code         ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "ThemisDustMix.hpp"
#include "DustEmGrainComposition.hpp"
#include "LogNormalGrainSizeDistribution.hpp"
#include "ModifiedPowerLawGrainSizeDistribution.hpp"

//////////////////////////////////////////////////////////////////////

void ThemisDustMix::setupSelfBefore()
{
    MultiGrainDustMix::setupSelfBefore();

    addPopulation(new DustEmGrainComposition(this, DustEmGrainComposition::GrainType::aOlM5, 2190.),
                  new LogNormalGrainSizeDistribution(this, 1.0e-9, 4900e-9, 8e-9, 1.), _numSilicateSizes,
                  GrainPopulation::NormalizationType::DustMassPerHydrogenMass, 0.255e-02);

    addPopulation(new DustEmGrainComposition(this, DustEmGrainComposition::GrainType::aPyM5, 2190.),
                  new LogNormalGrainSizeDistribution(this, 1.0e-9, 4900e-9, 8e-9, 1.), _numSilicateSizes,
                  GrainPopulation::NormalizationType::DustMassPerHydrogenMass, 0.255e-02);

    addPopulation(new DustEmGrainComposition(this, DustEmGrainComposition::GrainType::CM20, 1510.),
                  new LogNormalGrainSizeDistribution(this, 0.5e-9, 4900e-9, 7e-9, 1.), _numHydrocarbonSizes,
                  GrainPopulation::NormalizationType::DustMassPerHydrogenMass, 0.600e-03);

    addPopulation(new DustEmGrainComposition(this, DustEmGrainComposition::GrainType::CM20, 1600.),
                  new ModifiedPowerLawGrainSizeDistribution(this, 0.4e-9, 4900e-9, -5., 10e-9, 50e-9, 1., 1., 0., 1.),
                  _numHydrocarbonSizes, GrainPopulation::NormalizationType::DustMassPerHydrogenMass, 0.170e-02);
}

////////////////////////////////////////////////////////////////////
