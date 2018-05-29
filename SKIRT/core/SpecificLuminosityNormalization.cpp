/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "SpecificLuminosityNormalization.hpp"
#include "FatalError.hpp"
#include "SED.hpp"

////////////////////////////////////////////////////////////////////

double SpecificLuminosityNormalization::luminosity(SED* sed) const
{
    double Llambda = sed->specificLuminosity(_wavelength);
    if (Llambda <= 0) throw FATALERROR("The normalization wavelength is outside of the SED's wavelength range");
    return _specificLuminosity / Llambda;
}

//////////////////////////////////////////////////////////////////////
