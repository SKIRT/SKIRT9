/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef FILEBORDERWAVELENGTHGRID_HPP
#define FILEBORDERWAVELENGTHGRID_HPP

#include "DisjointWavelengthGrid.hpp"

////////////////////////////////////////////////////////////////////

/** FileBorderWavelengthGrid is a subclass of the DisjointWavelengthGrid class representing
    wavelength grids loaded from an input file. The floating point numbers in the first column of
    the text file specify the borders and possibly the characteristic wavelengths of the wavelength
    grid bins. The default unit is micron, but this can be overridden by a column header (see
    TextInFile). Any additional columns in the file are ignored.

    The precise behavior is governed by the value of the \em characteristic configuration property.
    If the value is \c Linear or \c Logarithmic, the wavelengths loaded from the file represent the
    borders of the wavelength grid bins. The order of the wavelengths in the file is not important;
    they will be sorted anyway. There must be at least two strictly positive wavelength values in
    the file, and duplicates are not allowed. The characteristic wavelength for each bin is derived
    automatically by linear or logarithmic interpolation depending on the value of \em
    characteristic.

    If the value of \em characteristic is \c Specified, the characteristic wavelengths for each bin
    are included in the wavelength list loaded from the file, interleaved with the border
    wavelengths (i.e., borders and characteristic wavelengths alternate). The number of values must
    be uneven and at least three. The list must be in strictly increasing or decreasing order,
    which means duplicates are not allowed, except that a zero characteristic wavelength indicates
    a segment that is not part of the grid, i.e. that lies between two non-adjacent bins. In other
    words, this option allows to (1) arbitrarily place characteristic wavelengths within each bin
    and (2) to specify intermediate wavelength ranges that are not covered by any bin. */
class FileBorderWavelengthGrid : public DisjointWavelengthGrid
{
    /** The enumeration type indicating how to determine the characteristic wavelength for each
        bin. See the class header for more information. */
    ENUM_DEF(Characteristic, Linear, Logarithmic, Specified)
        ENUM_VAL(Characteristic, Linear, "using linear interpolation")
        ENUM_VAL(Characteristic, Logarithmic, "using logarithmic interpolation")
        ENUM_VAL(Characteristic, Specified, "specified in the imported wavelength list")
    ENUM_END()

    ITEM_CONCRETE(FileBorderWavelengthGrid, DisjointWavelengthGrid,
                  "a wavelength grid loaded from a text file listing bin borders")
        ATTRIBUTE_TYPE_DISPLAYED_IF(FileBorderWavelengthGrid, "Level2")

        PROPERTY_STRING(filename, "the name of the file with the wavelength bin borders")

        PROPERTY_ENUM(characteristic, Characteristic, "determine the characteristic wavelength for each bin")
        ATTRIBUTE_DEFAULT_VALUE(characteristic, "Logarithmic")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function sets the wavelength grid to the bin borders loaded from the specified file.
        */
    void setupSelfBefore() override;
};

////////////////////////////////////////////////////////////////////

#endif
