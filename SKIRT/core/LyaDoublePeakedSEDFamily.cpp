/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "LyaDoublePeakedSEDFamily.hpp"
#include "Constants.hpp"
#include "NR.hpp"

//////////////////////////////////////////////////////////////////////

namespace
{
    // number of wavelength points in tabulated result per dispersion unit
    const int numWavelengthsPerDispersionUnit = 100;

    // double-peaked spectrum centered on 0 with scale of 1, evaluated at x
    double unitSpectrum(double x) { return 1.5 * x * x / (1. + cosh(x * x * x)); }
}

////////////////////////////////////////////////////////////////////

vector<SnapshotParameter> LyaDoublePeakedSEDFamily::parameterInfo() const
{
    return {SnapshotParameter::custom("line luminosity", "bolluminosity", "W"),
            SnapshotParameter::custom("scale", "velocity", "km/s")};
}

////////////////////////////////////////////////////////////////////

Range LyaDoublePeakedSEDFamily::intrinsicWavelengthRange() const
{
    return Range(1201e-10, 1231e-10);
}

////////////////////////////////////////////////////////////////////

double LyaDoublePeakedSEDFamily::specificLuminosity(double wavelength, const Array& parameters) const
{
    double L = parameters[0];
    double s = parameters[1];
    double wavelengthCenter = Constants::lambdaLya();
    double wavelengthDispersion = s * wavelengthCenter / Constants::c();
    return L * unitSpectrum((wavelength - wavelengthCenter) / wavelengthDispersion) / wavelengthDispersion;
}

////////////////////////////////////////////////////////////////////

double LyaDoublePeakedSEDFamily::cdf(Array& lambdav, Array& pv, Array& Pv, const Range& wavelengthRange,
                                     const Array& parameters) const
{
    double L = parameters[0];
    double s = parameters[1];
    double wavelengthCenter = Constants::lambdaLya();
    double wavelengthDispersion = s * wavelengthCenter / Constants::c();

    // build an appropriate grid
    size_t n = numWavelengthsPerDispersionUnit * wavelengthRange.width() / wavelengthDispersion;
    NR::buildLinearGrid(lambdav, wavelengthRange.min(), wavelengthRange.max(), n);

    // calculate the tabulated values
    pv.resize(n + 1);
    for (size_t i = 0; i <= n; ++i)
        pv[i] = L * unitSpectrum((lambdav[i] - wavelengthCenter) / wavelengthDispersion) / wavelengthDispersion;

    // calculate the cumulative distribution and normalization
    return NR::cdf2(false, lambdav, pv, Pv);
}

////////////////////////////////////////////////////////////////////
