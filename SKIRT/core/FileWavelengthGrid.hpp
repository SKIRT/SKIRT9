/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef FILEWAVELENGTHGRID_HPP
#define FILEWAVELENGTHGRID_HPP

#include "WavelengthGrid.hpp"

////////////////////////////////////////////////////////////////////

/** FileWavelengthGrid is a subclass of the WavelengthGrid class representing wavelength grids
    read  from a file. The text file must contain on the first line the number of wavelength grid
    points and subsequently all grid points in micron. The order of the grid points is irrelevant.
    */
class FileWavelengthGrid : public WavelengthGrid
{
    ITEM_CONCRETE(FileWavelengthGrid, WavelengthGrid, "a wavelength grid read from a file")

    PROPERTY_STRING(filename, "the name of the file with the wavelength grid points")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function reads the wavelength grid points from the specified file. */
    void setupSelfBefore() override;
};

////////////////////////////////////////////////////////////////////

#endif
