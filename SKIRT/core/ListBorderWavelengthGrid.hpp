/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef LISTBORDERWAVELENGTHGRID_HPP
#define LISTBORDERWAVELENGTHGRID_HPP

#include "DisjointWavelengthGrid.hpp"

//////////////////////////////////////////////////////////////////////

/** ListBorderWavelengthGrid is a subclass of the DisjointWavelengthGrid class that allows an
    arbitrary wavelength grid to be fully specified inside the configuration file (i.e. without
    referring to an input file). It is intended for use in cases where there are just a few
    wavelengths, but nothing keeps the user from specifying a long list.

    The precise behavior is governed by the value of the \em characteristic configuration property.
    If the value is \c Linear or \c Logarithmic, the wavelengths in the list represent the borders
    of the wavelength grid bins. The order of the wavelengths in the list is not important; they
    will be sorted anyway. There must be at least two strictly positive wavelength values in the
    list, and duplicates are not allowed. The characteristic wavelength for each bin is derived
    automatically by linear or logarithmic interpolation depending on the value of \em
    characteristic.

    If the value of \em characteristic is \c Specified, the characteristic wavelengths for each bin
    are included in the wavelength list, interleaved with the border wavelengths (i.e., borders and
    characteristic wavelengths alternate). The number of values must be uneven and at least three.
    The list must be in strictly increasing or decreasing order, which means duplicates are not
    allowed, except that a zero characteristic wavelength indicates a segment that is not part of
    the grid, i.e. that lies between two non-adjacent bins. In other words, this option allows to
    (1) arbitrarily place characteristic wavelengths within each bin and (2) to specify
    intermediate wavelength ranges that are not covered by any bin. */
class ListBorderWavelengthGrid : public DisjointWavelengthGrid
{
    /** The enumeration type indicating how to determine the characteristic wavelength for each
        bin. See the class header for more information. */
    ENUM_DEF(Characteristic, Linear, Logarithmic, Specified)
        ENUM_VAL(Characteristic, Linear, "using linear interpolation")
        ENUM_VAL(Characteristic, Logarithmic, "using logarithmic interpolation")
        ENUM_VAL(Characteristic, Specified, "specified in the configured wavelength list")
    ENUM_END()

    ITEM_CONCRETE(ListBorderWavelengthGrid, DisjointWavelengthGrid,
                  "a wavelength grid configured as a list of bin borders")
        ATTRIBUTE_TYPE_DISPLAYED_IF(ListBorderWavelengthGrid, "Level2")

        PROPERTY_DOUBLE_LIST(wavelengths, "the wavelength bin borders")
        ATTRIBUTE_QUANTITY(wavelengths, "wavelength")
        ATTRIBUTE_MIN_VALUE(wavelengths, "0 pm")
        ATTRIBUTE_MAX_VALUE(wavelengths, "1 m")

        PROPERTY_ENUM(characteristic, Characteristic, "determine the characteristic wavelength for each bin")
        ATTRIBUTE_DEFAULT_VALUE(characteristic, "Logarithmic")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function sets the wavelength grid to the bin borders specified by the \em wavelengths
        property. */
    void setupSelfBefore() override;
};

//////////////////////////////////////////////////////////////////////

#endif
