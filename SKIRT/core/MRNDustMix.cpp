/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "MRNDustMix.hpp"
#include "DraineGraphiteGrainComposition.hpp"
#include "DraineSilicateGrainComposition.hpp"
#include "PowerLawGrainSizeDistribution.hpp"

//////////////////////////////////////////////////////////////////////

void MRNDustMix::setupSelfBefore()
{
    MultiGrainDustMix::setupSelfBefore();

    // MRN grain size distribution coefficients taken from Weingartner & Draine (2001, ApJ, 548, 296) page 296
    const double amin = 5e-9;                  // 50 Angstrom
    const double amax = 250e-9;                // 0.25 micron
    const double gamma = 3.5;                  // power-law exponent
    const double Cs = pow(10, -25.11) * 1e-5;  // front factor for silicate converted from cm^2.5 to m^2.5
    const double Cg = pow(10, -25.13) * 1e-5;  // front factor for graphite converted from cm^2.5 to m^2.5

    auto sizeDistribution = new PowerLawGrainSizeDistribution(this, amin, amax, gamma);

    addPopulation(new DraineSilicateGrainComposition(this), sizeDistribution, _numSilicateSizes,
                  GrainPopulation::NormalizationType::FactorOnSizeDistribution, Cs);

    addPopulation(new DraineGraphiteGrainComposition(this), sizeDistribution, _numGraphiteSizes,
                  GrainPopulation::NormalizationType::FactorOnSizeDistribution, Cg);
}

//////////////////////////////////////////////////////////////////////
