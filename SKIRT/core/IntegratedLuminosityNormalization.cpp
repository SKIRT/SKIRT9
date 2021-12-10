/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "IntegratedLuminosityNormalization.hpp"
#include "FatalError.hpp"
#include "SED.hpp"

////////////////////////////////////////////////////////////////////

double IntegratedLuminosityNormalization::luminosityForSED(SED* sed) const
{
    // if the integration range is the source range, the configured luminosity is the answer
    if (_wavelengthRange == WavelengthRange::Source) return _integratedLuminosity;

    // determine the integration range (custom or all)
    double minWavelength = _wavelengthRange == WavelengthRange::Custom ? _minWavelength : 1e-10;
    double maxWavelength = _wavelengthRange == WavelengthRange::Custom ? _maxWavelength : 1;

    // refuse an emnpty wavelength range
    if (minWavelength >= maxWavelength) throw FATALERROR("the normalization wavelength range is empty");

    // get the normalized integrated luminosity from the sed
    double L = sed->integratedLuminosity(Range(minWavelength, maxWavelength));

    // refuse zero luminosity
    if (L <= 0) throw FATALERROR("the normalization luminosity is zero");

    return _integratedLuminosity / L;
}

//////////////////////////////////////////////////////////////////////
