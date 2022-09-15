/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "FileBorderWavelengthGrid.hpp"
#include "TextInFile.hpp"

////////////////////////////////////////////////////////////////////

void FileBorderWavelengthGrid::setupSelfBefore()
{
    DisjointWavelengthGrid::setupSelfBefore();

    // read the borders from the input file
    TextInFile infile(this, _filename, "wavelength grid");
    infile.addColumn("wavelength", "wavelength", "micron");
    Array wavelengths;
    infile.readAllColumns(wavelengths);
    infile.close();

    // set the wavelength grid
    switch (_characteristic)
    {
        case Characteristic::Linear: setWavelengthBorders(wavelengths, false); break;
        case Characteristic::Logarithmic: setWavelengthBorders(wavelengths, true); break;
        case Characteristic::Specified: setWavelengthSegments(wavelengths); break;
    }
}

//////////////////////////////////////////////////////////////////////
