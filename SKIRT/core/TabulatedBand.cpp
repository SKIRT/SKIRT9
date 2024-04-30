/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "TabulatedBand.hpp"
#include "FatalError.hpp"

////////////////////////////////////////////////////////////////////

void TabulatedBand::setupSelfBefore()
{
    Band::setupSelfBefore();

    // load the wavelengths and transmission values
    getWavelengthsAndTransmissions(_lambdav, _transv);

    // verify the number of values
    if (_lambdav.size() != _transv.size()) throw FATALERROR("Number of wavelength and transmission values differ");
    if (_lambdav.size() < 2) throw FATALERROR("Number of wavelength/transmission values is less than 2");

    // reverse the arrays if needed to get the wavelengths in increasing order
    if (_lambdav[0] > _lambdav[_lambdav.size() - 1])
    {
        std::reverse(begin(_lambdav), end(_lambdav));
        std::reverse(begin(_transv), end(_transv));
    }

    // verify that the wavelengths now are in increasing order
    if (!std::is_sorted(begin(_lambdav), end(_lambdav), std::less<>()))
        throw FATALERROR("Wavelengths are not sorted in increasing or decreasing order");

    // verify that the wavelength range is nonzero
    if (_lambdav[0] >= _lambdav[_lambdav.size() - 1]) throw FATALERROR("Wavelength range is empty");

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

size_t TabulatedBand::dataSize() const
{
    return _transv.size();
}

////////////////////////////////////////////////////////////////////

const double* TabulatedBand::wavelengthData() const
{
    return &_lambdav[0];
}

////////////////////////////////////////////////////////////////////

const double* TabulatedBand::transmissionData() const
{
    return &_transv[0];
}

////////////////////////////////////////////////////////////////////
