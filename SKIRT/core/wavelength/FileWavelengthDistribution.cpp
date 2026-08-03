/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "FileWavelengthDistribution.hpp"
#include "TextInFile.hpp"

//////////////////////////////////////////////////////////////////////

void FileWavelengthDistribution::getWavelengthsAndProbabilities(Array& lambdav, Array& pv) const
{
    // read the wavelengths and probabilities from the input file
    TextInFile infile(this, _filename, "wavelength probability distribution");
    infile.addColumn("wavelength", "wavelength", "micron");
    infile.addColumn("probability density", "specific", "W/m");
    infile.readAllColumns(lambdav, pv);
    infile.close();
}

//////////////////////////////////////////////////////////////////////
