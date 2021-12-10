/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef RESOLUTIONBORDERWAVELENGTHGRID_HPP
#define RESOLUTIONBORDERWAVELENGTHGRID_HPP

#include "DisjointWavelengthGrid.hpp"

////////////////////////////////////////////////////////////////////

/** ResolutionBorderWavelengthGrid is a subclass of the DisjointWavelengthGrid class representing
    logarithmically distributed wavelength grids with given spectral resolution (rather than number
    of bins) and given outer bin borders (rather than characteristic wavelengths).

    The number of bins is determined from the configured spectral resolution
    \f$R=\lambda/\Delta\lambda\f$ through \f[ N = \mathrm{ceil}\left[
    \frac{\log(\lambda_\mathrm{max}/\lambda_\mathrm{min})} {\log(1+1/R)} \right]. \f]

    The bin widths are equally distributed in log space between the specified minimum and maximum
    wavelength. The characteristic wavelength for each bin is determined as the geometric average
    of the bin's left and right borders. The grid always has at least one bin, which then has the
    specified minimum and maximum wavelength as its left and right border. */
class ResolutionBorderWavelengthGrid : public DisjointWavelengthGrid
{
    ITEM_CONCRETE(ResolutionBorderWavelengthGrid, DisjointWavelengthGrid,
                  "a logarithmic wavelength grid with given spectral resolution and outer borders")

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
    /** This function constructs the list of bin borders distributed logarithmically between
        \f$\lambda_{\text{min}}\f$ and \f$\lambda_{\text{max}}\f$. */
    void setupSelfBefore() override;
};

////////////////////////////////////////////////////////////////////

#endif
