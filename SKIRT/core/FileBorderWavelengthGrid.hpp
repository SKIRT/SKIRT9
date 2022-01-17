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
    the text file specify the borders of the wavelength grid bins. The default unit is micron, but
    this can be overridden by a column header (see TextInFile). Any additional columns in the file
    are ignored. The order of the wavelengths in the file is not important; they will be sorted
    anyway.

    The wavelengths loaded from the file represent the borders of the wavelength grid bins; there
    must be at least two borders in the file. The characteristic wavelength for each bin is derived
    automatically in linear or logarithmic space depending on the value of the \em log flag. For
    more details, refer to the DisjointWavelengthGrid::setWavelengthBorder() function. */
class FileBorderWavelengthGrid : public DisjointWavelengthGrid
{
    ITEM_CONCRETE(FileBorderWavelengthGrid, DisjointWavelengthGrid,
                  "a wavelength grid loaded from a text file listing bin borders")
        ATTRIBUTE_TYPE_DISPLAYED_IF(FileBorderWavelengthGrid, "Level2")

        PROPERTY_STRING(filename, "the name of the file with the wavelength bin borders")

        PROPERTY_BOOL(log, "use logarithmic scale")
        ATTRIBUTE_DEFAULT_VALUE(log, "true")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function sets the wavelength grid to the bin borders loaded from the specified file.
        */
    void setupSelfBefore() override;
};

////////////////////////////////////////////////////////////////////

#endif
