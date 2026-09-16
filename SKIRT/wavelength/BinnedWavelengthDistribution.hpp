/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef BINNEDWAVELENGTHDISTRIBUTION_HPP
#define BINNEDWAVELENGTHDISTRIBUTION_HPP

#include "GridWavelengthDistribution.hpp"

////////////////////////////////////////////////////////////////////

/** BinnedWavelengthDistribution is a specialty class for representing wavelength probability
    distributions derived from a wavelength grid that can be configured by the user, in which the
    generated wavelengths are spread across the full width of each bin instead of being
    concentrated at the bin's characteristic wavelength. This allows a relatively coarse wavelength
    grid to still yield a quasi-continuous sampling of wavelength space, for example to properly
    represent the continuum near an emission or absorption line without requiring an unnecessarily
    fine wavelength grid (and the corresponding increase in computational cost).

    When a BinnedWavelengthDistribution instance is used as the wavelength bias distribution for a
    source with a composite bias factor of one, the source will emit photon packets across the full
    width of the (in-range) bins of the configured grid, with equal probability per bin, and with a
    uniform probability density within each bin. More precisely, the BinnedWavelengthDistribution
    instance uses only the bins of the configured grid for which the characteristic wavelength falls
    inside the wavelength range of the associated source (obtained through the
    SourceWavelengthRangeInterface). If none of the characteristic wavelengths fall inside the
    source range, a fatal error is issued.

    See also DiscreteWavelengthDistribution, which offers similar functionality but concentrates
    the generated wavelengths at the bin's characteristic wavelength instead. This is required, for
    example, when mimicking codes that emit photon packets at a fixed set of discrete wavelengths. */
class BinnedWavelengthDistribution : public GridWavelengthDistribution
{
    ITEM_CONCRETE(BinnedWavelengthDistribution, GridWavelengthDistribution,
                  "a binned wavelength probability distribution derived from a wavelength grid")
    ITEM_END()

    //======================== Other Functions =======================

protected:
    /** This function returns a random wavelength, uniformly distributed across the full width of
        the given bin. */
    double wavelengthInBin(int ell) const override;
};

////////////////////////////////////////////////////////////////////

#endif
