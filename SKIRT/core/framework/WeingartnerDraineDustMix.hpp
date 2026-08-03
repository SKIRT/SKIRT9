/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef WEINGARTNERDRAINEDUSTMIX_HPP
#define WEINGARTNERDRAINEDUSTMIX_HPP

#include "MultiGrainDustMix.hpp"

////////////////////////////////////////////////////////////////////

/** The WeingartnerDraineDustMix class represents realistic dust mixtures consisting of populations
    of silicate, graphite, and PAH dust grains. It is assumed that 50\% of the PAHs are neutral and
    50\% are ionized. The size distributions of each of these grains are fitted in such a way that
    they can reproduce the extinction curve of the Milky Way, the LMC, or the SMC, respectively.
    For details we refer to Weingartner & Draine (2001, ApJ, 548, 296) and Li & Draine (2001, ApJ,
    554, 778).

    As noted at the end of section 3 in Weingartner & Draine (2001), the absence of PAHs in the SMC
    dust model might be a special attribute of the specific sight line, so it may be unrealistic.

    For each of the populations, the user can configure the number of grain size bins used to
    discretize the thermal emission calculations for the dust mix.

    The class uses the grain compositions defined in the DraineSilicateGrainComposition,
    DraineGraphiteGrainComposition, DraineNeutralPAHGrainComposition, and
    DraineIonizedPAHGrainComposition classes.

    The grain size distributions are implementated internally in the class because they are very
    specific to this particular dust mix. They are given as complicated analytical formulas with
    predefined constant parameters depending on the selected environment (Milky Way, LMC or SMC).
    For the graphite and silicate populations, the distribution is a power-law function with a
    curvature and an exponential cutoff. For the PAH populations, it is the sum of two log-normal
    distributions, cut off at some upper grain size. The exact forms of the grain size distribution
    functions can be found in Weingartner & Draine (2001, ApJ, 548, 296). */
class WeingartnerDraineDustMix : public MultiGrainDustMix
{
    /** The enumeration type indicating the environment for the Weingartner-Draine dust. */
    ENUM_DEF(Environment, MilkyWay, LMC, SMC)
        ENUM_VAL(Environment, MilkyWay, "the Milky Way")
        ENUM_VAL(Environment, LMC, "the Large Magellanic Cloud")
        ENUM_VAL(Environment, SMC, "the Small Magellanic Cloud")
    ENUM_END()

    ITEM_CONCRETE(WeingartnerDraineDustMix, MultiGrainDustMix, "a Weingartner and Draine (2001) dust mix")

        PROPERTY_ENUM(environment, Environment, "the environment determining the dust model")
        ATTRIBUTE_DEFAULT_VALUE(environment, "MilkyWay")

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
