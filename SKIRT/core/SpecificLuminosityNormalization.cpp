/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "SpecificLuminosityNormalization.hpp"
#include "ContSED.hpp"
#include "FatalError.hpp"
#include "Units.hpp"

////////////////////////////////////////////////////////////////////

double SpecificLuminosityNormalization::luminosityForSED(SED* sed) const
{
    auto contsed = dynamic_cast<ContSED*>(sed);
    if (!contsed) throw FATALERROR("Cannot use specific luminosity normalization for a line emission spectrum");

    // get the normalized specific luminosity in the SED
    double LlambdaSED = contsed->specificLuminosity(_wavelength);
    if (LlambdaSED <= 0) throw FATALERROR("The normalization wavelength is outside of the SED's wavelength range");

    // convert the user-configured specific luminosity to per-wavelength units
    double LlambdaUser = Units::fromFluxStyle(_wavelength, _specificLuminosity, Units::fluxStyle(_unitStyle));

    // return the ratio
    return LlambdaUser / LlambdaSED;
}

//////////////////////////////////////////////////////////////////////
