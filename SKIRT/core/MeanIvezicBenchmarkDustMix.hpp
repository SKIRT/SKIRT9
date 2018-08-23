/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef MEANIVEZICBENCHMARKDUSTMIX_HPP
#define MEANIVEZICBENCHMARKDUSTMIX_HPP

#include "SingleGrainDustMix.hpp"

////////////////////////////////////////////////////////////////////

/** The MeanIvezicBenchmarkDustMix class represents the idealized dust mix used in the 1D radiative
    transfer benchmark calculations of Ivezic et al. (1997, MNRAS, 291, 121). Scattering is assumed
    to be isotropic and the absorption and scattering coefficients are approximated by simple
    analytical functions.

    With \f$\lambda\f$ representing the wavelength expressed in \f$\mu\f$m, the absorption
    coefficient is given by \f[ \frac{\kappa_\lambda^{\text{abs}}} {\kappa_1^{\text{abs}}} =
    \begin{cases} \; 1 & \text{if $\lambda<1$} \\ \; \dfrac{1}{\lambda} & \text{else}, \end{cases}
    \f] and the scattering coefficient by \f[ \frac{\kappa_\lambda^{\text{sca}}}
    {\kappa_1^{\text{sca}}} = \begin{cases} \; 1 & \text{if $\lambda<1$} \\ \; \dfrac{1}{\lambda^4}
    & \text{else}. \end{cases} \f]

    The extinction coefficients in the benchmark data are scale-free; we arbitrarily scale them to
    a reasonable order of magnitude. */
class MeanIvezicBenchmarkDustMix : public SingleGrainDustMix
{
    ITEM_CONCRETE(MeanIvezicBenchmarkDustMix, SingleGrainDustMix,
                  "a dust mix from the Ivezic 1D benchmark (mean properties)")
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
