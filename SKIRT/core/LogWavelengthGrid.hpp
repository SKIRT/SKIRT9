/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef LOGWAVELENGTHGRID_HPP
#define LOGWAVELENGTHGRID_HPP

#include "DisjointWavelengthGrid.hpp"

////////////////////////////////////////////////////////////////////

/** LogWavelengthGrid is a subclass of the DisjointWavelengthGrid class representing
    logarithmically distributed wavelength grids. The characteristic wavelengths of the grid bins
    are equally distributed (in log space) between and including the specified minimum and maximum
    wavelength. The outermost bins are given the same width as the inner bins (in log space), which
    implies that the outermost bin borders are placed beyond the specified minimum and maximum
    wavelength. The grid must have at least two bins, which then have the specified minimum and
    maximum wavelength as their respective characteristic wavelength. */
class LogWavelengthGrid : public DisjointWavelengthGrid
{
    ITEM_CONCRETE(LogWavelengthGrid, DisjointWavelengthGrid, "a logarithmic wavelength grid")

        PROPERTY_DOUBLE(minWavelength, "the shortest wavelength")
        ATTRIBUTE_QUANTITY(minWavelength, "wavelength")
        ATTRIBUTE_MIN_VALUE(minWavelength, "1 pm")
        ATTRIBUTE_MAX_VALUE(minWavelength, "1 m")

        PROPERTY_DOUBLE(maxWavelength, "the longest wavelength")
        ATTRIBUTE_QUANTITY(maxWavelength, "wavelength")
        ATTRIBUTE_MIN_VALUE(maxWavelength, "1 pm")
        ATTRIBUTE_MAX_VALUE(maxWavelength, "1 m")

        PROPERTY_INT(numWavelengths, "the number of wavelength grid points")
        ATTRIBUTE_MIN_VALUE(numWavelengths, "2")
        ATTRIBUTE_DEFAULT_VALUE(numWavelengths, "25")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function constructs the list of \f$N\f$ characteristic wavelengths distributed
        logarithmically between \f$\lambda_{\text{min}}\f$ and \f$\lambda_{\text{max}}\f$. */
    void setupSelfBefore() override;
};

////////////////////////////////////////////////////////////////////

#endif
