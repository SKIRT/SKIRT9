/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "FileBand.hpp"
#include "TextInFile.hpp"

////////////////////////////////////////////////////////////////////

void FileBand::getWavelengthsAndTransmissions(Array& lambdav, Array& transv) const
{
    // read the wavelengths and transmission values from the input file
    TextInFile infile(this, _filename, "transmission curve");
    infile.addColumn("wavelength", "wavelength", "micron");
    infile.addColumn("transmission value");
    infile.readAllColumns(lambdav, transv);
    infile.close();
}

////////////////////////////////////////////////////////////////////
