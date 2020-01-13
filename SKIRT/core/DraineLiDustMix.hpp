/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef DRAINELIDUSTMIX_HPP
#define DRAINELIDUSTMIX_HPP

#include "MultiGrainDustMix.hpp"

////////////////////////////////////////////////////////////////////

/** The DraineLiDustMix class represents a dust mixture of silicate, graphite, and PAH dust grains
    designed by Draine \& Li 2007. The size distribution of each of these dust grain populations is
    fine-tuned in such a way that the global dust properties accurately reproduce the extinction
    curve of the Milky Way.

    For each of the populations, the user can configure the number of grain size bins used to
    discretize the thermal emission calculations for the dust mix.

    The following is an extract from the GRAIN_DL07.DAT data file specifying this dust model for the
    DustEM code (version 4.2), described by Compiègne et al. 2011 (AA, 525, A103):

    \verbatim
    # grain type, nsize, type keywords, Mdust/MH, rho, amin, amax, alpha/a0 [, at, ac, gamma (ED)] [, au, zeta, eta (CV)]
    # (cgs units)
    aSil      70 plaw-ed-cv 7.64E-03 3.50E+0 3.10E-8 2.00E-4 -3.21E+0 1.64E-5 1.00E-5 3E+0 1.64E-5  3.00E-01  1E+0
    Gra       70 plaw-ed-cv 2.21E-03 2.24E+0 3.10E-8 2.00E-4 -2.54E+0 1.07E-6 4.28E-5 3E+0 1.07E-6 -1.65E-01  1E+0
    Gra       30 logn       1.66E-04 2.24E+0 3.10E-8 4.00E-6  2.00E-7 5.50E-1
    PAH0_DL07 10 mix-logn   4.97E-04 2.24E+0 3.10E-8 1.20E-7  4.00E-8 4.00E-1
    PAH1_DL07 10 mix-logn   4.97E-04 2.24E+0 3.10E-8 1.20E-7  4.00E-8 4.00E-1
    \endverbatim

    This class uses grain compositions obtained through the DustEmGrainComposition class,
    respectively selecting the grain types aSil (astronimical silicates), Gra (graphite), PAH0_DL07
    (neutral PAHs), and PAH1_DL07 (ionized PAHs). It further uses the size distributions defined in
    the ModifiedPowerLawGrainSizeDistribution and LogNormalGrainSizeDistribution classes. */
class DraineLiDustMix : public MultiGrainDustMix
{
    ITEM_CONCRETE(DraineLiDustMix, MultiGrainDustMix, "a Draine and Li (2007) dust mix")

        PROPERTY_INT(numSilicateSizes, "the number of silicate grain size bins")
        ATTRIBUTE_MIN_VALUE(numSilicateSizes, "1")
        ATTRIBUTE_DEFAULT_VALUE(numSilicateSizes, "5")

        PROPERTY_INT(numGraphiteSizes, "the number of graphite grain size bins (for each of two populations)")
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
