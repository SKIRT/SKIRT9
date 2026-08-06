/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef FILEWAVELENGTHGRID_HPP
#define FILEWAVELENGTHGRID_HPP

#include "DisjointWavelengthGrid.hpp"

////////////////////////////////////////////////////////////////////

/** FileWavelengthGrid is a subclass of the DisjointWavelengthGrid class representing wavelength
    grids loaded from an input file. The floating point numbers in the first column of the text
    file specify the characteristic wavelengths. The default unit is micron, but this can be
    overridden by a column header (see TextInFile). Any additional columns in the file are ignored.
    The order of the wavelengths in the file is not important; they will be sorted anyway.

    In all cases, the wavelengths loaded from the file represent the characteristic wavelengths of
    the grid bins, and the bin borders are derived automatically. Note that the outermost bin
    borders are placed beyond the outermost characteristic wavelengths.

    If \em relativeHalfWidth is zero (the default value) the class constructs a consecutive range
    of adjacent wavelength bins in linear or logarithmic space depending on the value of the \em
    log flag. For more details, refer to the DisjointWavelengthGrid::setWavelengthRange() function.

    If \em relativeHalfWidth is nonzero, the class constructs a set of distinct nonadjacent
    wavelength bins with the specified relative half bin width. In this case, the value of the \em
    log flag is ignored. For more details, refer to the DisjointWavelengthGrid::setWavelengthBins()
    function. */
class FileWavelengthGrid : public DisjointWavelengthGrid
{
    ITEM_CONCRETE(FileWavelengthGrid, DisjointWavelengthGrid, "a wavelength grid loaded from a text file")
        ATTRIBUTE_TYPE_DISPLAYED_IF(FileWavelengthGrid, "Level2")

        PROPERTY_STRING(filename, "the name of the file with the characteristic wavelengths")

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
    /** This function reads the wavelength grid points from the specified file. */
    void setupSelfBefore() override;
};

////////////////////////////////////////////////////////////////////

#endif
