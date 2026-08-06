/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef ZUBKODUSTMIX_HPP
#define ZUBKODUSTMIX_HPP

#include "MultiGrainDustMix.hpp"

////////////////////////////////////////////////////////////////////

/** The ZubkoDustMix class represents a realistic dust mixture of bare (i.e. non-composite)
    silicate, graphite, neutral PAH and ionized PAH dust grains. The size distribution of each of
    these dust grain populations is finetuned in such a way that the global dust properties
    accurately reproduce the extinction, emission and abundance constraints on the Milky Way. The
    size distributions are taken from Zubko, Dwek & Arendt (2004, ApJS, 152, 211), model BARE_GR_S.
    It is assumed that 50% of the PAH grains are neutral and 50% are ionized.

    For each of the populations, the user can configure the number of grain size bins used to
    discretize the thermal emission calculations for the dust mix.

    The class uses the grain compositions defined in the DraineSilicateGrainComposition,
    DraineGraphiteGrainComposition, DraineNeutralPAHGrainComposition, and
    DraineIonizedPAHGrainComposition classes, and the size distributions defined in
    the ZubkoSilicateGrainSizeDistribution, ZubkoGraphiteGrainSizeDistribution, and
    ZubkoPAHGrainSizeDistribution classes. */
class ZubkoDustMix : public MultiGrainDustMix
{
    ITEM_CONCRETE(ZubkoDustMix, MultiGrainDustMix, "a Zubko et al. (2004) dust mix")

        PROPERTY_INT(numSilicateSizes, "the number of silicate grain size bins")
        ATTRIBUTE_MIN_VALUE(numSilicateSizes, "1")
        ATTRIBUTE_DEFAULT_VALUE(numSilicateSizes, "5")

        PROPERTY_INT(numGraphiteSizes, "the number of graphite grain size bins")
        ATTRIBUTE_MIN_VALUE(numGraphiteSizes, "1")
        ATTRIBUTE_DEFAULT_VALUE(numGraphiteSizes, "5")

        PROPERTY_INT(numPAHSizes, "the number of neutral and ionized PAH size bins (each)")
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
