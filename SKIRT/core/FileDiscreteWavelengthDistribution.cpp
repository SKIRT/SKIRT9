/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "FileDiscreteWavelengthDistribution.hpp"
#include "TextInFile.hpp"

//////////////////////////////////////////////////////////////////////

Array FileDiscreteWavelengthDistribution::getWavelengths() const
{
    // read the wavelengths from the input file
    Array lambdav;
    TextInFile infile(this, _filename, "wavelengths");
    infile.addColumn("wavelength", "wavelength", "micron");
    infile.readAllColumns(lambdav);
    infile.close();

    return lambdav;
}

//////////////////////////////////////////////////////////////////////
