/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef MEANINTERSTELLARDUSTMIX_HPP
#define MEANINTERSTELLARDUSTMIX_HPP

#include "SingleGrainDustMix.hpp"

////////////////////////////////////////////////////////////////////

/** The MeanInterstellarDustMix class represents a population of identical dust grains with
    properties approximating those of a mixture that is appropriate for the typical interstellar
    dust medium. The model consists of a mixture of carbonaceous grains and amorphous silicate
    grains. Carbonaceous grains are PAH-like when small, and graphite-like when large. (see Li &
    Draine 2001, ApJ, 554, 778). Size distributions are taken from Weingartner & Draine (2001, ApJ,
    548, 296), more specifically, they refer to the case A model for \f$R_V=3.1\f$, renormalized
    following Draine (2003, ApJ, 598, 1017). Specifically, the grain abundances relative to H have
    been reduced by a factor 0.93 relative to the \f$R_V=3.1\f$ size distribution with
    \f$\text{[C/H]}_{\text{PAH}} = 60~{\text{ppm}}\f$ in Weingartner & Draine (2001, ApJ, 548,
    296). The PAH C abundance relative to H is assumed to be \f$\text{[C/H]}_{\text{PAH}} = 0.93
    \times 60~{\text{ppm}} = 55.8~{\text{ppm}}\f$.

    The dust mass per hydrogen nucleon is obtained from information in the header of the data file.
    The optical data were downloaded from Bruce Draine's home page
    https://www.astro.princeton.edu/~draine/ and more specifically from
    ftp://ftp.astro.princeton.edu/draine/dust/mix/kext_albedo_WD_MW_3.1_60_D03.all */
class MeanInterstellarDustMix : public SingleGrainDustMix
{
    ITEM_CONCRETE(MeanInterstellarDustMix, SingleGrainDustMix, "a typical interstellar dust mix (mean properties)")
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
