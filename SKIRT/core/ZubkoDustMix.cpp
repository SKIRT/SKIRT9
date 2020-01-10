/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "ZubkoDustMix.hpp"
#include "DraineGraphiteGrainComposition.hpp"
#include "DraineIonizedPAHGrainComposition.hpp"
#include "DraineNeutralPAHGrainComposition.hpp"
#include "DraineSilicateGrainComposition.hpp"
#include "GrainPopulation.hpp"
#include "ZubkoGraphiteGrainSizeDistribution.hpp"
#include "ZubkoPAHGrainSizeDistribution.hpp"
#include "ZubkoSilicateGrainSizeDistribution.hpp"

////////////////////////////////////////////////////////////////////

void ZubkoDustMix::setupSelfBefore()
{
    MultiGrainDustMix::setupSelfBefore();

    addPopulation(new DraineSilicateGrainComposition(this), new ZubkoSilicateGrainSizeDistribution(this),
                  _numSilicateSizes, GrainPopulation::NormalizationType::FactorOnSizeDistribution, 1.);

    addPopulation(new DraineGraphiteGrainComposition(this), new ZubkoGraphiteGrainSizeDistribution(this),
                  _numGraphiteSizes, GrainPopulation::NormalizationType::FactorOnSizeDistribution, 1.);

    addPopulation(new DraineNeutralPAHGrainComposition(this), new ZubkoPAHGrainSizeDistribution(this), _numPAHSizes,
                  GrainPopulation::NormalizationType::FactorOnSizeDistribution, 0.5);

    addPopulation(new DraineIonizedPAHGrainComposition(this), new ZubkoPAHGrainSizeDistribution(this), _numPAHSizes,
                  GrainPopulation::NormalizationType::FactorOnSizeDistribution, 0.5);
}

////////////////////////////////////////////////////////////////////
