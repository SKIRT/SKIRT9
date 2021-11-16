/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "SpecificLuminosityNormalization.hpp"
#include "ContSED.hpp"
#include "FatalError.hpp"
#include "Units.hpp"

////////////////////////////////////////////////////////////////////

double SpecificLuminosityNormalization::luminosity(SED* sed) const
{
    auto contsed = dynamic_cast<ContSED*>(sed);
    if (!contsed) throw FATALERROR("Cannot use specific luminosity normalization for a line emission spectrum");

    // get the normalized specific luminosity in the SED
    double LlambdaSED = contsed->specificLuminosity(_wavelength);
    if (LlambdaSED <= 0) throw FATALERROR("The normalization wavelength is outside of the SED's wavelength range");

    // convert the user-configured specific luminosity to per-wavelength units
    double LlambdaUser = 0.;
    switch (_unitStyle)
    {
        case UnitStyle::neutralmonluminosity:
            LlambdaUser = Units::fromNeutralStyle(_wavelength, _specificLuminosity);
            break;
        case UnitStyle::wavelengthmonluminosity:
            LlambdaUser = Units::fromWavelengthStyle(_wavelength, _specificLuminosity);
            break;
        case UnitStyle::frequencymonluminosity:
            LlambdaUser = Units::fromFrequencyStyle(_wavelength, _specificLuminosity);
            break;
        case UnitStyle::energymonluminosity:
            LlambdaUser = Units::fromEnergyStyle(_wavelength, _specificLuminosity);
            break;
    }

    // return the ratio
    return LlambdaUser / LlambdaSED;
}

//////////////////////////////////////////////////////////////////////
