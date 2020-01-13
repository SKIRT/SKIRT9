/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef TRUSTBENCHMARKDUSTMIX_HPP
#define TRUSTBENCHMARKDUSTMIX_HPP

#include "MultiGrainDustMix.hpp"

////////////////////////////////////////////////////////////////////

/** The TrustBenchmarkDustMix class and represents a realistic dust mixture of bare (i.e.
    non-composite) silicate, graphite, and PAH dust grains. The size distribution of each of these
    dust grain populations is fine-tuned in such a way that the global dust properties accurately
    reproduce the extinction, emission and abundance constraints on the Milky Way. The size
    distributions are taken from Zubko, Dwek & Arendt (2004, ApJS, 152, 211) and correspond to
    model BARE_GR_S.

    For each of the populations, the user can configure the number of grain size bins used to
    discretize the thermal emission calculations for the dust mix.

    This dust mix has been prepared by Karl Misselt for the TRUST benchmark simulations; see Camps
    et al. 2015 (A&A 580:A87) and Gordon et al. 2017 (A&A 603:A114). The data can be downloaded
    from http://www.shg.ugent.be/html/_downloads.html or http://ipag.osug.fr/RT13/RTTRUST/opa.php

    The class uses the grain compositions defined in the TrustSilicateGrainComposition,
    TrustGraphiteGrainComposition, and TrustNeutralPAHGrainComposition, and the size distributions
    defined in the ZubkoSilicateGrainSizeDistribution, ZubkoGraphiteGrainSizeDistribution, and
    ZubkoPAHGrainSizeDistribution classes. */
class TrustBenchmarkDustMix : public MultiGrainDustMix
{
    ITEM_CONCRETE(TrustBenchmarkDustMix, MultiGrainDustMix, "a TRUST benchmark dust mix")
        ATTRIBUTE_TYPE_DISPLAYED_IF(TrustBenchmarkDustMix, "Level2")

        PROPERTY_INT(numSilicateSizes, "the number of silicate grain size bins")
        ATTRIBUTE_MIN_VALUE(numSilicateSizes, "1")
        ATTRIBUTE_DEFAULT_VALUE(numSilicateSizes, "5")

        PROPERTY_INT(numGraphiteSizes, "the number of graphite grain size bins")
        ATTRIBUTE_MIN_VALUE(numGraphiteSizes, "1")
        ATTRIBUTE_DEFAULT_VALUE(numGraphiteSizes, "5")

        PROPERTY_INT(numPAHSizes, "the number of PAH size bins")
        ATTRIBUTE_MIN_VALUE(numPAHSizes, "1")
        ATTRIBUTE_DEFAULT_VALUE(numPAHSizes, "5")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function adds the relevant grain populations to the dust mix */
    void setupSelfBefore() override;
};

////////////////////////////////////////////////////////////////////

#endif
