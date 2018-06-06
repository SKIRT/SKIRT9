/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "FileWavelengthGrid.hpp"
#include "TextInFile.hpp"

////////////////////////////////////////////////////////////////////

void FileWavelengthGrid::setupSelfBefore()
{
    WavelengthGrid::setupSelfBefore();

    // read the wavelengths from the input file
    TextInFile infile(this, _filename, "wavelength grid");
    const vector<Array>& columns = infile.readAllColumns(1);
    infile.close();

    // set the wavelength grid (convert from micron to m)
    if (_relativeHalfWidth) setWavelengthBins(columns[0]*1e-6, _relativeHalfWidth);
    else setWavelengthRange(columns[0]*1e-6);
}

//////////////////////////////////////////////////////////////////////
