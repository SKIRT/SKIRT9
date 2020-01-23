/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "LyaGaussianSED.hpp"
#include "Constants.hpp"
#include "FatalError.hpp"
#include "NR.hpp"
#include "Random.hpp"

//////////////////////////////////////////////////////////////////////

namespace
{
    // intrinsic half-range in dispersion units
    const double intrinsicRange = 9;

    // number of wavelength points in tabulated result per dispersion unit
    const int numWavelengthsPerDispersionUnit = 100;

    // Gaussian centered on 0 with dispersion of 1, evaluated at x
    double unitGaussian(double x) { return (0.5 * M_SQRT1_2 * M_2_SQRTPI) * exp(-0.5 * x * x); }

    // integrate unit Gaussian from x1 to x2
    double integrateUnitGaussian(double x1, double x2) { return 0.5 * (erf(M_SQRT1_2 * x2) - erf(M_SQRT1_2 * x1)); }
}

//////////////////////////////////////////////////////////////////////

void LyaGaussianSED::setupSelfBefore()
{
    SED::setupSelfBefore();

    _wavelengthCenter = Constants::lambdaLya();
    _wavelengthDispersion = _wavelengthCenter * dispersion() / Constants::c();
    _wavelengthRange.set(_wavelengthCenter - intrinsicRange * _wavelengthDispersion,
                         _wavelengthCenter + intrinsicRange * _wavelengthDispersion);

    // verify that the source wavelength range contains the intrinsic range so we can normalize to the intrinsic range
    if (!normalizationWavelengthRange().contains(_wavelengthRange))
        throw FATALERROR("The source wavelength range must fully contain the intrinsic range of this SED");
}

//////////////////////////////////////////////////////////////////////

Range LyaGaussianSED::intrinsicWavelengthRange() const
{
    return _wavelengthRange;
}

//////////////////////////////////////////////////////////////////////

double LyaGaussianSED::specificLuminosity(double wavelength) const
{
    return _wavelengthRange.contains(wavelength)
               ? unitGaussian((wavelength - _wavelengthCenter) / _wavelengthDispersion) / _wavelengthDispersion
               : 0.;
}

//////////////////////////////////////////////////////////////////////

void LyaGaussianSED::specificLuminosityArray(Array& lambdav, Array& pv, const Range& wavelengthRange) const
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
        size_t n = max(static_cast<size_t>(10), static_cast<size_t>(numWavelengthsPerDispersionUnit
                                                                    * intersection.width() / _wavelengthDispersion));
        NR::buildLinearGrid(lambdav, intersection.min(), intersection.max(), n);

        // calculate the tabulated values
        pv.resize(n + 1);
        for (size_t i = 0; i <= n; ++i)
            pv[i] = unitGaussian((lambdav[i] - _wavelengthCenter) / _wavelengthDispersion) / _wavelengthDispersion;
    }
}

//////////////////////////////////////////////////////////////////////

double LyaGaussianSED::integratedLuminosity(const Range& wavelengthRange) const
{
    // if the given range includes the complete intrinsic range, the result is trivial
    if (wavelengthRange.contains(_wavelengthRange)) return 1.;

    // calculate the intersection between the given range and our intrinsic range
    Range intersection = _wavelengthRange;
    intersection.intersect(wavelengthRange);

    // if the intersection is empty, the result is trivial
    if (intersection.empty()) return 0.;

    // otherwise, integrate over the intersection
    return integrateUnitGaussian((intersection.min() - _wavelengthCenter) / _wavelengthDispersion,
                                 (intersection.max() - _wavelengthCenter) / _wavelengthDispersion);
}

//////////////////////////////////////////////////////////////////////

double LyaGaussianSED::generateWavelength() const
{
    while (true)
    {
        double wavelength = random()->gauss() * _wavelengthDispersion + _wavelengthCenter;
        if (_wavelengthRange.contains(wavelength)) return wavelength;
    }
}

//////////////////////////////////////////////////////////////////////
