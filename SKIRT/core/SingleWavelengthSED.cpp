/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "SingleWavelengthSED.hpp"
#include "FatalError.hpp"

//////////////////////////////////////////////////////////////////////

Range SingleWavelengthSED::intrinsicWavelengthRange() const
{
    return Range(_wavelength * (1 - 1e-6), _wavelength * (1 + 1e-6));  // avoid an empty range
}

//////////////////////////////////////////////////////////////////////

double SingleWavelengthSED::specificLuminosity(double /*wavelength*/) const
{
    throw FATALERROR("Specific luminosity is undefined for SingleWavelengthSED");
}

//////////////////////////////////////////////////////////////////////

void SingleWavelengthSED::specificLuminosityArray(Array& /*lambdav*/, Array& /*pv*/,
                                                  const Range& /*wavelengthRange*/) const
{
    throw FATALERROR("Specific luminosity array is undefined for SingleWavelengthSED");
}

//////////////////////////////////////////////////////////////////////

double SingleWavelengthSED::integratedLuminosity(const Range& wavelengthRange) const
{
    return wavelengthRange.contains(_wavelength) ? 1. : 0.;
}

//////////////////////////////////////////////////////////////////////

double SingleWavelengthSED::generateWavelength() const
{
    return _wavelength;
}

//////////////////////////////////////////////////////////////////////
