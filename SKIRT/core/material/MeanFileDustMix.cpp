/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "MeanFileDustMix.hpp"
#include "TextInFile.hpp"

//////////////////////////////////////////////////////////////////////

double MeanFileDustMix::getDustProperties(Array& lambdav, Array& kappaextv, Array& albedov, Array& asymmparv) const
{
    // read the wavelengths and optical properties from the input file
    TextInFile infile(this, _filename, "optical dust properties");
    infile.addColumn("wavelength", "wavelength", "micron");
    infile.addColumn("extinction mass coefficient", "masscoefficient", "m2/kg");
    infile.addColumn("scattering albedo");
    infile.addColumn("scattering asymmetry parameter");
    infile.readAllColumns(lambdav, kappaextv, albedov, asymmparv);
    infile.close();

    return 1.5e-29;  // in kg/H
}

//////////////////////////////////////////////////////////////////////
