/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "FileBand.hpp"
#include "FatalError.hpp"
#include "TextInFile.hpp"

////////////////////////////////////////////////////////////////////

void FileBand::setupSelfBefore()
{
    Band::setupSelfBefore();

    // read the wavelengths and transmission values from the input file
    TextInFile infile(this, _filename, "transmission curve");
    infile.addColumn("wavelength", "wavelength", "micron");
    infile.addColumn("transmission value");
    infile.readAllColumns(_lambdav, _transv);
    infile.close();

    // verify the number of values
    if (_transv.size() < 2) throw FATALERROR("Number of loaded transmission values is less than 2");

    // calculate the normalization factor for the transmission values
    size_t size = _transv.size();
    double norm = 0.;
    for (size_t i = 1; i != size; ++i)
    {
        double dlambda = _lambdav[i] - _lambdav[i - 1];
        double trans = 0.5 * (_transv[i - 1] + _transv[i]);
        norm += trans * dlambda;
    }

    // normalize the transmission values
    _transv /= norm;
}

////////////////////////////////////////////////////////////////////

size_t FileBand::dataSize() const
{
    return _transv.size();
}

////////////////////////////////////////////////////////////////////

const double* FileBand::wavelengthData() const
{
    return &_lambdav[0];
}

////////////////////////////////////////////////////////////////////

const double* FileBand::transmissionData() const
{
    return &_transv[0];
}

////////////////////////////////////////////////////////////////////
