/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "DraineLiDustMix.hpp"
#include "DustEmGrainComposition.hpp"
#include "LogNormalGrainSizeDistribution.hpp"
#include "ModifiedPowerLawGrainSizeDistribution.hpp"

////////////////////////////////////////////////////////////////////

void DraineLiDustMix::setupSelfBefore()
{
    MultiGrainDustMix::setupSelfBefore();

    addPopulation(
        new DustEmGrainComposition(this, DustEmGrainComposition::GrainType::aSil, 3500.),
        new ModifiedPowerLawGrainSizeDistribution(this, 3.1e-10, 2.0e-6, -3.21, 1.64e-7, 1.00e-7, 3., 1.64e-7, 0.3, 1.),
        _numSilicateSizes, GrainPopulation::NormalizationType::DustMassPerHydrogenMass, 7.64e-3);

    addPopulation(new DustEmGrainComposition(this, DustEmGrainComposition::GrainType::Gra, 2240.),
                  new ModifiedPowerLawGrainSizeDistribution(this, 3.1e-10, 2.0e-6, -2.54, 1.07e-8, 4.28e-7, 3., 1.07e-8,
                                                            -0.165, 1.),
                  _numGraphiteSizes, GrainPopulation::NormalizationType::DustMassPerHydrogenMass, 2.21e-3);

    addPopulation(new DustEmGrainComposition(this, DustEmGrainComposition::GrainType::Gra, 2240.),
                  new LogNormalGrainSizeDistribution(this, 3.1e-10, 4.0e-8, 2.0e-9, 0.55), _numGraphiteSizes,
                  GrainPopulation::NormalizationType::DustMassPerHydrogenMass, 1.66e-4);

    addPopulation(new DustEmGrainComposition(this, DustEmGrainComposition::GrainType::PAH0DL07, 2240.),
                  new LogNormalGrainSizeDistribution(this, 3.1e-10, 1.2e-9, 4.0e-10, 0.4), _numPAHSizes,
                  GrainPopulation::NormalizationType::DustMassPerHydrogenMass, 0.5 * 4.97e-4);

    addPopulation(new DustEmGrainComposition(this, DustEmGrainComposition::GrainType::PAH1DL07, 2240.),
                  new LogNormalGrainSizeDistribution(this, 3.1e-10, 1.2e-9, 4.0e-10, 0.4), _numPAHSizes,
                  GrainPopulation::NormalizationType::DustMassPerHydrogenMass, 0.5 * 4.97e-4);
}

////////////////////////////////////////////////////////////////////
