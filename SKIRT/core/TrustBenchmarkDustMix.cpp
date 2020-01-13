/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "TrustBenchmarkDustMix.hpp"
#include "TrustGraphiteGrainComposition.hpp"
#include "TrustNeutralPAHGrainComposition.hpp"
#include "TrustSilicateGrainComposition.hpp"
#include "ZubkoGraphiteGrainSizeDistribution.hpp"
#include "ZubkoPAHGrainSizeDistribution.hpp"
#include "ZubkoSilicateGrainSizeDistribution.hpp"

//////////////////////////////////////////////////////////////////////

void TrustBenchmarkDustMix::setupSelfBefore()
{
    MultiGrainDustMix::setupSelfBefore();

    addPopulation(new TrustSilicateGrainComposition(this), new ZubkoSilicateGrainSizeDistribution(this),
                  _numSilicateSizes, GrainPopulation::NormalizationType::FactorOnSizeDistribution, 1.);

    addPopulation(new TrustGraphiteGrainComposition(this), new ZubkoGraphiteGrainSizeDistribution(this),
                  _numGraphiteSizes, GrainPopulation::NormalizationType::FactorOnSizeDistribution, 1.);

    addPopulation(new TrustNeutralPAHGrainComposition(this), new ZubkoPAHGrainSizeDistribution(this), _numPAHSizes,
                  GrainPopulation::NormalizationType::FactorOnSizeDistribution, 1.);
}

////////////////////////////////////////////////////////////////////
