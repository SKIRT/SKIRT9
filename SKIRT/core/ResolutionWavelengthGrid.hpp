/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef RESOLUTIONWAVELENGTHGRID_HPP
#define RESOLUTIONWAVELENGTHGRID_HPP

#include "DisjointWavelengthGrid.hpp"

////////////////////////////////////////////////////////////////////

/** ResolutionWavelengthGrid is a subclass of the DisjointWavelengthGrid class representing
    logarithmically distributed wavelength grids with given spectral resolution (rather than number
    of bins) and given minimum and maximum characteristic wavelengths.

    Because the outer bin borders extend beyond the specified minimum and maximum wavelengths,
    there is an extra bin compared to the ResolutionBorderWavelengthGrid. In other words, the
    number of bins is determined from the configured spectral resolution
    \f$R=\lambda/\Delta\lambda\f$ through \f[ N = 1+\mathrm{ceil}\left[
    \frac{\log(\lambda_\mathrm{max}/\lambda_\mathrm{min})} {\log(1+1/R)} \right]. \f]

    The characteristic wavelengths of the grid bins are equally distributed (in log space) between
    and including the specified minimum and maximum wavelength. The outermost bins are given the
    same width as the inner bins (in log space), which implies that the outermost bin borders are
    placed beyond the specified minimum and maximum wavelength. The grid has at least two bins,
    which then have the specified minimum and maximum wavelength as their respective characteristic
    wavelength. */
class ResolutionWavelengthGrid : public DisjointWavelengthGrid
{
    ITEM_CONCRETE(ResolutionWavelengthGrid, DisjointWavelengthGrid,
                  "a logarithmic wavelength grid with given spectral resolution")

        PROPERTY_DOUBLE(minWavelength, "the shortest wavelength")
        ATTRIBUTE_QUANTITY(minWavelength, "wavelength")
        ATTRIBUTE_MIN_VALUE(minWavelength, "1 pm")
        ATTRIBUTE_MAX_VALUE(minWavelength, "1 m")

        PROPERTY_DOUBLE(maxWavelength, "the longest wavelength")
        ATTRIBUTE_QUANTITY(maxWavelength, "wavelength")
        ATTRIBUTE_MIN_VALUE(maxWavelength, "1 pm")
        ATTRIBUTE_MAX_VALUE(maxWavelength, "1 m")

        PROPERTY_DOUBLE(resolution, "the spectral resolution R of the grid")
        ATTRIBUTE_MIN_VALUE(resolution, "[1e-6")
        ATTRIBUTE_MAX_VALUE(resolution, "1e6]")
        ATTRIBUTE_DEFAULT_VALUE(resolution, "10")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function constructs the list of characteristic wavelengths distributed logarithmically
        between \f$\lambda_{\text{min}}\f$ and \f$\lambda_{\text{max}}\f$. */
    void setupSelfBefore() override;
};

////////////////////////////////////////////////////////////////////

#endif
