/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "LyaDoublePeakedSED.hpp"
#include "Constants.hpp"
#include "FatalError.hpp"
#include "NR.hpp"
#include "Random.hpp"

//////////////////////////////////////////////////////////////////////

namespace
{
    // intrinsic half-range in scale units
    const double intrinsicRange = 3.6;

    // number of wavelength points in tabulated result per scale unit
    const int numWavelengthsPerScaleUnit = 100;

    // double-peaked spectrum centered on 0 with scale of 1, evaluated at x
    double unitSpectrum(double x) { return 1.5 * x * x / (1. + cosh(x * x * x)); }

    // integrate double-peaked unit spectrum from x1 to x2
    double integrateUnitSpectrum(double x1, double x2)
    {
        return (1. / (1. + exp(-x2 * x2 * x2))) - (1. / (1. + exp(-x1 * x1 * x1)));
    }
}

//////////////////////////////////////////////////////////////////////

void LyaDoublePeakedSED::setupSelfBefore()
{
    SED::setupSelfBefore();

    _wavelengthCenter = Constants::lambdaLya();
    _wavelengthScale = _wavelengthCenter * scale() / Constants::c();
    _wavelengthRange.set(_wavelengthCenter - intrinsicRange * _wavelengthScale,
                         _wavelengthCenter + intrinsicRange * _wavelengthScale);

    // verify that the source wavelength range contains the intrinsic range so we can normalize to the intrinsic range
    if (!normalizationWavelengthRange().contains(_wavelengthRange))
        throw FATALERROR("The source wavelength range must fully contain the intrinsic range of this SED");
}

//////////////////////////////////////////////////////////////////////

Range LyaDoublePeakedSED::intrinsicWavelengthRange() const
{
    return _wavelengthRange;
}

//////////////////////////////////////////////////////////////////////

double LyaDoublePeakedSED::specificLuminosity(double wavelength) const
{
    return _wavelengthRange.contains(wavelength)
               ? unitSpectrum((wavelength - _wavelengthCenter) / _wavelengthScale) / _wavelengthScale
               : 0.;
}

//////////////////////////////////////////////////////////////////////

void LyaDoublePeakedSED::specificLuminosityArray(Array& lambdav, Array& pv, const Range& wavelengthRange) const
{
    // calculate the intersection between the given range and our intrinsic range
    Range intersection = _wavelengthRange;
    intersection.intersect(wavelengthRange);

    // if the intersection is empty, return an empty result
    if (intersection.empty())
    {
        lambdav.resize(0);
        pv.resize(0);
    }
    // otherwise, return a table with appropriate resolution
    else
    {
        // build an appropriate grid
        size_t n = max(static_cast<size_t>(10),
                       static_cast<size_t>(numWavelengthsPerScaleUnit * intersection.width() / _wavelengthScale));
        NR::buildLinearGrid(lambdav, intersection.min(), intersection.max(), n);

        // calculate the tabulated values
        pv.resize(n + 1);
        for (size_t i = 0; i <= n; ++i)
            pv[i] = unitSpectrum((lambdav[i] - _wavelengthCenter) / _wavelengthScale) / _wavelengthScale;
    }
}

//////////////////////////////////////////////////////////////////////

double LyaDoublePeakedSED::integratedLuminosity(const Range& wavelengthRange) const
{
    // if the given range includes the complete intrinsic range, the result is trivial
    if (wavelengthRange.contains(_wavelengthRange)) return 1.;

    // calculate the intersection between the given range and our intrinsic range
    Range intersection = _wavelengthRange;
    intersection.intersect(wavelengthRange);

    // if the intersection is empty, the result is trivial
    if (intersection.empty()) return 0.;

    // otherwise, integrate over the intersection
    return integrateUnitSpectrum((intersection.min() - _wavelengthCenter) / _wavelengthScale,
                                 (intersection.max() - _wavelengthCenter) / _wavelengthScale);
}

//////////////////////////////////////////////////////////////////////

double LyaDoublePeakedSED::generateWavelength() const
{
    while (true)
    {
        double uniform = random()->uniform();
        double sample = cbrt(log(uniform / (1 - uniform)));
        double wavelength = sample * _wavelengthScale + _wavelengthCenter;
        if (_wavelengthRange.contains(wavelength)) return wavelength;
    }
}

//////////////////////////////////////////////////////////////////////
