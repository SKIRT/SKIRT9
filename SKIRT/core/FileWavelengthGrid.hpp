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
    file specify the characteristic wavelengths in micron. Any additional columns in the file are
    ignored. The order of the wavelengths in the file is not important; they will be sorted anyway.

    The class constructs either a consecutive range of adjacent wavelength bins (when \em
    relativeHalfWidth is zero, the default value) or a set of distinct nonadjacent wavelength bins
    (with the relative half bin width given by a nonzero value of \em relativeHalfWidth). Refer to
    the DisjointWavelengthGrid class for more details. */
class FileWavelengthGrid : public DisjointWavelengthGrid
{
    ITEM_CONCRETE(FileWavelengthGrid, DisjointWavelengthGrid, "a wavelength grid loaded from a text file")
        ATTRIBUTE_TYPE_DISPLAYED_IF(FileWavelengthGrid, "Level2")

        PROPERTY_STRING(filename, "the name of the file with the characteristic wavelengths")

        PROPERTY_DOUBLE(relativeHalfWidth, "the relative half width for discrete bins, or zero for a consecutive range")
        ATTRIBUTE_MIN_VALUE(relativeHalfWidth, "[0")
        ATTRIBUTE_MAX_VALUE(relativeHalfWidth, "1[")
        ATTRIBUTE_DEFAULT_VALUE(relativeHalfWidth, "0")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function reads the wavelength grid points from the specified file. */
    void setupSelfBefore() override;
};

////////////////////////////////////////////////////////////////////

#endif
