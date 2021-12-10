/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef LINBORDERWAVELENGTHGRID_HPP
#define LINBORDERWAVELENGTHGRID_HPP

#include "DisjointWavelengthGrid.hpp"

////////////////////////////////////////////////////////////////////

/** LinBorderWavelengthGrid is a subclass of the DisjointWavelengthGrid class representing linearly
    distributed wavelength grids with given outer bin borders (rather than characteristic
    wavelengths).

    The bin widths are equally distributed between the specified minimum and maximum wavelength.
    The characteristic wavelength for each bin is determined as the arithmetic average of the bin's
    left and right borders. The grid must have at least one bin, which then has the specified
    minimum and maximum wavelength as its left and right border. */
class LinBorderWavelengthGrid : public DisjointWavelengthGrid
{
    ITEM_CONCRETE(LinBorderWavelengthGrid, DisjointWavelengthGrid, "a linear wavelength grid with given outer borders")

        PROPERTY_DOUBLE(minWavelength, "the shortest wavelength")
        ATTRIBUTE_QUANTITY(minWavelength, "wavelength")
        ATTRIBUTE_MIN_VALUE(minWavelength, "1 pm")
        ATTRIBUTE_MAX_VALUE(minWavelength, "1 m")

        PROPERTY_DOUBLE(maxWavelength, "the longest wavelength")
        ATTRIBUTE_QUANTITY(maxWavelength, "wavelength")
        ATTRIBUTE_MIN_VALUE(maxWavelength, "1 pm")
        ATTRIBUTE_MAX_VALUE(maxWavelength, "1 m")

        PROPERTY_INT(numWavelengthBins, "the number of wavelength bins")
        ATTRIBUTE_MIN_VALUE(numWavelengthBins, "1")
        ATTRIBUTE_DEFAULT_VALUE(numWavelengthBins, "25")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function constructs the list of \f$N+1\f$ bin borders distributed linearly between
        \f$\lambda_{\text{min}}\f$ and \f$\lambda_{\text{max}}\f$. */
    void setupSelfBefore() override;
};

////////////////////////////////////////////////////////////////////

#endif
