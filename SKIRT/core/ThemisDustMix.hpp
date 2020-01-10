/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef THEMISDUSTMIX_HPP
#define THEMISDUSTMIX_HPP

#include "MultiGrainDustMix.hpp"

////////////////////////////////////////////////////////////////////

/** The ThemisDustMix class represents the THEMIS model for dust in the diffuse interstellar medium
    described by Jones et al. 2017 (A&A, 602, A46) and the references therein. In this model, there
    are two families of dust particles: amorphous silicate and amorphous hydrocarbon. For the
    silicate, it is assumed that 50\% of the mass is amorphous enstatite, and that the remaining
    half is amorphous forsterite, and the size distribution is a lognormal distribution (the same
    for both populations). For the amorphous hydrocarbon, the size distribution is a combination of
    a lognormal and a power-law distribution.

    For each of the populations, the user can configure the number of grain size bins used to
    discretize the thermal emission calculations for the dust mix.

    The following is an extract from the GRAIN_J17.DAT file specifying the THEMIS dust model for
    the DustEM code, described by Compiègne et al. 2011 (AA, 525, A103):

    \verbatim
    # grain type, nsize, type keywords, Mdust/MH, rho, amin, amax, alpha/a0 [, at, ac, gamma (ED)] [, au, zeta, eta (CV)]
    # cgs units
    aOLM5  100 logn      0.255E-02  2.190E+00  1.00E-07   4900.0E-07     8.00E-07    1.00E+00
    aPyM5  100 logn      0.255E-02  2.190E+00  1.00E-07   4900.0E-07     8.00E-07    1.00E+00
    CM20   100 logn      0.600E-03  1.510E+00  0.50E-07   4900.0E-07     7.00E-07    1.00E+00
    CM20   100 plaw-ed   0.170E-02  1.600E+00  0.40E-07   4900.0E-07    -5.00E-00   10.00E-07   50.0E-07  1.000E+00
    \endverbatim

    The aOlM5 (olivine) grain type represents the amorphous silicates with forsterite-normative
    composition, aPyM5 (pyroxene) represents the amorphous silicates with enstatite-normative
    composition, and CM20 represents the amorphous carbonaceous dust grains in the THEMIS model.

    This class uses grain compositions obtained through the DustEmGrainComposition class,
    respectively selecting the grain types aOLM5, aPyM5, and CM20. It further uses the size
    distributions defined in the ModifiedPowerLawGrainSizeDistribution and
    LogNormalGrainSizeDistribution classes. */
class ThemisDustMix : public MultiGrainDustMix
{
    ITEM_CONCRETE(ThemisDustMix, MultiGrainDustMix, "a THEMIS (Jones et al. 2017) dust mix")

        PROPERTY_INT(numSilicateSizes, "the number of grain size bins for each of the silicate populations")
        ATTRIBUTE_MIN_VALUE(numSilicateSizes, "1")
        ATTRIBUTE_DEFAULT_VALUE(numSilicateSizes, "5")

        PROPERTY_INT(numHydrocarbonSizes, "the number of grain size bins for each of the hydrocarbon populations")
        ATTRIBUTE_MIN_VALUE(numHydrocarbonSizes, "1")
        ATTRIBUTE_DEFAULT_VALUE(numHydrocarbonSizes, "5")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function adds the relevant grain populations to the dust mix */
    void setupSelfBefore() override;
};

////////////////////////////////////////////////////////////////////

#endif
