/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "LineLuminosityNormalization.hpp"
#include "FatalError.hpp"
#include "SED.hpp"

////////////////////////////////////////////////////////////////////

double LineLuminosityNormalization::luminosityForSED(SED* sed) const
{
    // get the normalized luminosity for the configured line from the sed, using a very narrow wavelength range
    double L = sed->integratedLuminosity(Range(_wavelength * (1 - 1e-10), _wavelength * (1 + 1e-10)));

    // refuse zero luminosity
    if (L <= 0) throw FATALERROR("the configured normalization wavelength does not correspond to an emission line");

    return _luminosity / L;
}

//////////////////////////////////////////////////////////////////////
