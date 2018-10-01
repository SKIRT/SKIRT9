/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "TabulatedSED.hpp"
#include "Configuration.hpp"
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

    // determine the intersected wavelength range
    Range range(_inlambdav[0], _inlambdav[_inlambdav.size()-1]);
    range.intersect(find<Configuration>()->sourceWavelengthRange());
    if (range.empty()) throw FATALERROR("SED wavelength range does not overlap source wavelength range");

    // construct the regular and cumulative distributions in the intersected range
    double norm = NR::cdf<NR::interpolateLogLog>(_lambdav, _pv, _Pv, _inlambdav, _inpv, range);

    // also normalize the intrinsic distribution
    _inpv /= norm;
}

//////////////////////////////////////////////////////////////////////

double TabulatedSED::specificLuminosity(double wavelength) const
{
    int i = NR::locateFail(_inlambdav, wavelength);
    if (i < 0) return 0.;
    return NR::interpolateLogLog(wavelength, _inlambdav[i], _inlambdav[i+1], _inpv[i], _inpv[i+1]);
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
