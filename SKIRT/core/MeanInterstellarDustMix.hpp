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
    ftp://ftp.astro.princeton.edu/draine/dust/mix/kext_albedo_WD_MW_3.1_60_D03.all

    <b>Extreme forward scattering for X-ray wavelengths</b>

    The properties for this dust mix as discussed above are given down to wavelengths as short as
    \f$10^{-4}\f$ micron, i.e. well into the soft X-ray wavelength range. The tabulated cross
    sections seem to be very accurate, however the average scattering angle cosine values in this
    regime are unrealistic (e.g., equal to one for the shortest wavelengths).

    Draine 2003c (ApJ, 98, 1026–1037) presents experimental data for the phase function in the
    soft X-ray wavelength range and proposes an analytical phase function to model the observed
    extreme forward scattering. Rather than attempting to implement this phase function, we
    approximate it through the Henyey-Greenstein phase function just as we do for other
    wavelengths. For wavelengths shorter than 0.01 micron, the value of the HG asymmetry parameter
    \f$g\f$ that causes the HG function to best approximate the Draine 2003c function for extreme
    forward angles can be written as \f$g = 1 - \lambda/(\mu\mathrm{m})\f$, i.e. one minus the
    wavelength in micron. We thus replace the asymmetry parameter values in that wavelength range
    by this improved appoximation before storing the result in the SKIRT resource table. */
class MeanInterstellarDustMix : public SingleGrainDustMix
{
    ITEM_CONCRETE(MeanInterstellarDustMix, SingleGrainDustMix, "a typical interstellar dust mix (mean properties)")
    ITEM_END()

    //======================== Other Functions =======================

protected:
    /** This function returns the name of the stored table resource tabulating the basic optical
        properties for this dust mix. */
    string resourceNameForOpticalProps() const override;
};

////////////////////////////////////////////////////////////////////

#endif
