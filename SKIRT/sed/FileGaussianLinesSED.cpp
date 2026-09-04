/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "FileGaussianLinesSED.hpp"
#include "TextInFile.hpp"

//////////////////////////////////////////////////////////////////////

void FileGaussianLinesSED::getLineProperties(Array& lambdav, Array& dispersionv, Array& Lv) const
{
    // read the wavelengths, velocity dispersions and luminosities from the input file
    TextInFile infile(this, _filename, "Gaussian line definitions");
    infile.addColumn("wavelength", "wavelength", "micron");
    infile.addColumn("velocity dispersion", "velocity", "km/s");
    infile.addColumn("luminosity", "bolluminosity", "Lsun");
    infile.readAllColumns(lambdav, dispersionv, Lv);
    infile.close();
}

//////////////////////////////////////////////////////////////////////
