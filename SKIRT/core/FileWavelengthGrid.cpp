/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "FileWavelengthGrid.hpp"
#include "TextInFile.hpp"

////////////////////////////////////////////////////////////////////

void FileWavelengthGrid::setupSelfBefore()
{
    DisjointWavelengthGrid::setupSelfBefore();

    // read the wavelengths from the input file
    TextInFile infile(this, _filename, "wavelength grid");
    infile.addColumn("wavelength", "wavelength", "micron");
    Array wavelengths;
    infile.readAllColumns(wavelengths);
    infile.close();

    // set the wavelength grid
    if (_relativeHalfWidth)
        setWavelengthBins(wavelengths, _relativeHalfWidth);
    else
        setWavelengthRange(wavelengths, _log);
}

//////////////////////////////////////////////////////////////////////
