/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef MEANDRAINELIDUSTMIX_HPP
#define MEANDRAINELIDUSTMIX_HPP

#include "SingleGrainDustMix.hpp"

////////////////////////////////////////////////////////////////////

/** The MeanDraineLiDustMix class represents a population of identical dust grains with properties
    approximating those of a mixture of graphite, silicate and PAH dust grains designed by Draine
    \& Li 2007. The size distribution of each of these dust grain populations is fine-tuned in such
    a way that the global dust properties accurately reproduce the extinction curve of the Milky
    Way. The optical properties of this dust mixture are generated using the DustEm code using the
    following input data:

    \verbatim
    # grain type, nsize, type keywords, Mdust/MH, rho, amin, amax, alpha/a0 [, at, ac, gamma (ED)] [, au, zeta, eta (CV)]
    # (cgs units)
    PAH0_DL07 10 mix-logn   5.40E-4 2.24E+0 3.10E-08 1.20E-7  3.50E-8 4.00E-1
    PAH1_DL07 10 mix-logn   5.40E-4 2.24E+0 3.10E-08 1.20E-7  3.50E-8 4.00E-1
    Gra       30 logn       1.80E-4 2.24E+0 3.10E-08 4.00E-6  3.00E-7 4.00E-1
    Gra       70 plaw-ed-cv 2.33E-3 2.24E+0 3.10E-08 2.00E-4 -2.54E+0 1.07E-6 4.28E-5 3.00E+0 1.07E-6 -1.65E-1 1.00E+0
    aSil      70 plaw-ed-cv 8.27E-3 3.50E+0 3.10E-08 2.00E-4 -3.21E+0 1.64E-5 1.00E-5 3.00E+0 1.64E-5  3.00E-1 1.00E+0
    \endverbatim
*/
class MeanDraineLiDustMix : public SingleGrainDustMix
{
    ITEM_CONCRETE(MeanDraineLiDustMix, SingleGrainDustMix,
                  "a Draine & Li 2007 dust mix (mean properties)")
    ITEM_END()

    //======================== Other Functions =======================

    /** This function returns the scattering mode supported by this material mix. For this dust
        mix, it always returns the basic HenyeyGreenstein scattering mode. */
    ScatteringMode scatteringMode() const override;

    /** This function returns the name of the stored table resource tabulating the basic optical
        properties for this dust mix. */
    string resourceNameForOpticalProps() const override;
};

////////////////////////////////////////////////////////////////////

#endif
