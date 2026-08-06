/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef LISTWAVELENGTHGRID_HPP
#define LISTWAVELENGTHGRID_HPP

#include "DisjointWavelengthGrid.hpp"

//////////////////////////////////////////////////////////////////////

/** ListWavelengthGrid is a subclass of the DisjointWavelengthGrid class that allows an arbitrary
    wavelength grid to be fully specified inside the configuration file (i.e. without referring to
    an input file). It is intended for use in cases where there are just a few wavelengths, but
    nothing keeps the user from specifying a long list. The order of the specified wavelengths is
    not important; they will be sorted anyway.

    In all cases, the wavelengths specified in the list represent the characteristic wavelengths of
    the grid bins, and the bin borders are derived automatically. Note that the outermost bin
    borders are placed beyond the outermost characteristic wavelengths.

    If \em relativeHalfWidth is zero (the default value) the class constructs a consecutive range
    of adjacent wavelength bins in linear or logarithmic space depending on the value of the \em
    log flag. For more details, refer to the DisjointWavelengthGrid::setWavelengthRange() function.

    If \em relativeHalfWidth is nonzero, the class constructs a set of distinct nonadjacent
    wavelength bins with the specified relative half bin width. In this case, the value of the \em
    log flag is ignored. For more details, refer to the DisjointWavelengthGrid::setWavelengthBins()
    function. */
class ListWavelengthGrid : public DisjointWavelengthGrid
{
    ITEM_CONCRETE(ListWavelengthGrid, DisjointWavelengthGrid, "a wavelength grid configured as a list")
        ATTRIBUTE_TYPE_DISPLAYED_IF(ListWavelengthGrid, "Level2")

        PROPERTY_DOUBLE_LIST(wavelengths, "the characteristic wavelength for each bin")
        ATTRIBUTE_QUANTITY(wavelengths, "wavelength")
        ATTRIBUTE_MIN_VALUE(wavelengths, "1 pm")
        ATTRIBUTE_MAX_VALUE(wavelengths, "1 m")

        PROPERTY_DOUBLE(relativeHalfWidth, "the relative half width for discrete bins, or zero for a consecutive range")
        ATTRIBUTE_MIN_VALUE(relativeHalfWidth, "[0")
        ATTRIBUTE_MAX_VALUE(relativeHalfWidth, "1[")
        ATTRIBUTE_DEFAULT_VALUE(relativeHalfWidth, "0")

        PROPERTY_BOOL(log, "use logarithmic scale")
        ATTRIBUTE_DEFAULT_VALUE(log, "true")
        ATTRIBUTE_RELEVANT_IF(log, "!relativeHalfWidth")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function sets the wavelength grid to the characteristic wavelengths specified by the
        \em wavelengths property. */
    void setupSelfBefore() override;
};

//////////////////////////////////////////////////////////////////////

#endif
