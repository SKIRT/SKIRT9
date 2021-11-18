/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "TabulatedSED.hpp"
#include "FatalError.hpp"
#include "NR.hpp"
#include "Random.hpp"

//////////////////////////////////////////////////////////////////////

void TabulatedSED::setupSelfBefore()
{
    SED::setupSelfBefore();

    // obtain the wavelengths and specific luminosities
    getWavelengthsAndLuminosities(_inlambdav, _inpv);

    // verify number of points
    if (_inlambdav.size() < 2) throw FATALERROR("SED must have at least two wavelength/luminosity pairs");

    // reverse the arrays if needed to get the wavelengths in increasing order
    if (_inlambdav[0] > _inlambdav[_inlambdav.size() - 1])
    {
        std::reverse(begin(_inlambdav), end(_inlambdav));
        std::reverse(begin(_inpv), end(_inpv));
    }

    // construct the regular and cumulative distributions
    double norm = NR::cdf<NR::interpolateLogLog>(_lambdav, _pv, _Pv, _inlambdav, _inpv, normalizationWavelengthRange());

    // also normalize the intrinsic distribution
    _inpv /= norm;
}

//////////////////////////////////////////////////////////////////////

Range TabulatedSED::intrinsicWavelengthRange() const
{
    return Range(_inlambdav[0], _inlambdav[_inlambdav.size() - 1]);
}

//////////////////////////////////////////////////////////////////////

double TabulatedSED::specificLuminosity(double wavelength) const
{
    return NR::value<NR::interpolateLogLog>(wavelength, _inlambdav, _inpv);
}

//////////////////////////////////////////////////////////////////////

void TabulatedSED::specificLuminosityArray(Array& lambdav, Array& pv, const Range& wavelengthRange) const
{
    Array Pv;  // the contents of this array is not used, so this could be optimized if needed
    NR::cdf<NR::interpolateLogLog>(lambdav, pv, Pv, _inlambdav, _inpv, wavelengthRange);
}

//////////////////////////////////////////////////////////////////////

double TabulatedSED::integratedLuminosity(const Range& wavelengthRange) const
{
    Array lambdav, pv, Pv;  // the contents of these arrays is not used, so this could be optimized if needed
    return NR::cdf<NR::interpolateLogLog>(lambdav, pv, Pv, _inlambdav, _inpv, wavelengthRange);
}

//////////////////////////////////////////////////////////////////////

double TabulatedSED::generateWavelength() const
{
    return random()->cdfLogLog(_lambdav, _pv, _Pv);
}

//////////////////////////////////////////////////////////////////////
