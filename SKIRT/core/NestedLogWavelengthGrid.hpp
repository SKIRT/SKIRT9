/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef NESTEDLOGWAVELENGTHGRID_HPP
#define NESTEDLOGWAVELENGTHGRID_HPP

#include "DisjointWavelengthGrid.hpp"

////////////////////////////////////////////////////////////////////

/** NestedLogWavelengthGrid is a subclass of the DisjointWavelengthGrid class representing hybrid
    grids consisting of a logarithmically distributed wavelength grid in which another, more
    compact logarithmic grid is embedded. It can be very useful to get higher-resolution spectra in
    a particular wavelength grid while still covering a broad wavelength range.

    The characteristic wavelengths of the bins in each of the two grids (low-resolution grid and
    high-resolution subgrid) are equally distributed (in log space) between and including the
    specified minimum and maximum wavelength. The outermost bins are given the same width as the
    inner bins (in log space), which implies that the outermost bin borders are placed beyond the
    specified minimum and maximum wavelength. The grids must have at least two bins, which then
    have the specified minimum and maximum wavelength as their respective characteristic
    wavelength. */
class NestedLogWavelengthGrid : public DisjointWavelengthGrid
{
    ITEM_CONCRETE(NestedLogWavelengthGrid, DisjointWavelengthGrid, "a nested logarithmic wavelength grid")

        PROPERTY_DOUBLE(minWavelengthBaseGrid, "the shortest wavelength of the low-resolution grid")
        ATTRIBUTE_QUANTITY(minWavelengthBaseGrid, "wavelength")
        ATTRIBUTE_MIN_VALUE(minWavelengthBaseGrid, "1 pm")
        ATTRIBUTE_MAX_VALUE(minWavelengthBaseGrid, "1 m")

        PROPERTY_DOUBLE(maxWavelengthBaseGrid, "the longest wavelength of the low-resolution grid")
        ATTRIBUTE_QUANTITY(maxWavelengthBaseGrid, "wavelength")
        ATTRIBUTE_MIN_VALUE(maxWavelengthBaseGrid, "1 pm")
        ATTRIBUTE_MAX_VALUE(maxWavelengthBaseGrid, "1 m")

        PROPERTY_INT(numWavelengthsBaseGrid, "the number of wavelength grid points in the low-resolution grid")
        ATTRIBUTE_MIN_VALUE(numWavelengthsBaseGrid, "2")
        ATTRIBUTE_DEFAULT_VALUE(numWavelengthsBaseGrid, "25")

        PROPERTY_DOUBLE(minWavelengthSubGrid, "the shortest wavelength of the high-resolution subgrid")
        ATTRIBUTE_QUANTITY(minWavelengthSubGrid, "wavelength")
        ATTRIBUTE_MIN_VALUE(minWavelengthSubGrid, "1 pm")
        ATTRIBUTE_MAX_VALUE(minWavelengthSubGrid, "1 m")

        PROPERTY_DOUBLE(maxWavelengthSubGrid, "the longest wavelength of the high-resolution subgrid")
        ATTRIBUTE_QUANTITY(maxWavelengthSubGrid, "wavelength")
        ATTRIBUTE_MIN_VALUE(maxWavelengthSubGrid, "1 pm")
        ATTRIBUTE_MAX_VALUE(maxWavelengthSubGrid, "1 m")

        PROPERTY_INT(numWavelengthsSubGrid, "the number of wavelength grid points in the high-resolution subgrid")
        ATTRIBUTE_MIN_VALUE(numWavelengthsSubGrid, "2")
        ATTRIBUTE_DEFAULT_VALUE(numWavelengthsSubGrid, "25")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function constructs the list of characteristic wavelengths. First, \f$N\f$ wavelength
        grid points are distributed logarithmically between \f$\lambda_{\text{min}}\f$ and
        \f$\lambda_{\text{max}}\f$. Next, \f$N_{\text{zoom}}\f$ wavelength grid points are
        distributed logarithmically between \f$\lambda_{\text{zoom,min}}\f$ and
        \f$\lambda_{\text{zoom,max}}\f$. Both sets of wavelength grids are subsequently merged and
        overlapping grid points are removed. */
    void setupSelfBefore() override;
};

////////////////////////////////////////////////////////////////////

#endif
