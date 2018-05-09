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
    const vector<Array>& rows = infile.readAllRows(1);
    infile.close();

    // copy the result, converting from micron to m
    size_t n = rows.size();
    Array lambdav(n);
    for (size_t i=0; i<n; i++) lambdav[i] = rows[i][0] * 1e-6;

    // set the wavelength grid
    if (_relativeHalfWidth) setWavelengthBins(lambdav, _relativeHalfWidth);
    else setWavelengthRange(lambdav);
}

//////////////////////////////////////////////////////////////////////
