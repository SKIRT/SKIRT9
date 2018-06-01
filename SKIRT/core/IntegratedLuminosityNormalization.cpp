/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "IntegratedLuminosityNormalization.hpp"
#include "SED.hpp"

////////////////////////////////////////////////////////////////////

double IntegratedLuminosityNormalization::luminosity(SED* sed) const
{
    // if the integration range is the source range, the configured luminosity is the answer
    if (_wavelengthRange == WavelengthRange::Source) return _integratedLuminosity;

    // determine the integration range (custom or all)
    double minWavelength = _wavelengthRange == WavelengthRange::Custom ? _minWavelength : 1e-10;
    double maxWavelength = _wavelengthRange == WavelengthRange::Custom ? _maxWavelength : 1;

    // get the normalized integrated luminosity from the sed
    double L = sed->integratedLuminosity(minWavelength, maxWavelength);
    return _integratedLuminosity / L;
}

//////////////////////////////////////////////////////////////////////
