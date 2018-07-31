/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "BandLuminosityNormalization.hpp"
#include "Constants.hpp"
#include "FatalError.hpp"
#include "SED.hpp"

////////////////////////////////////////////////////////////////////

double BandLuminosityNormalization::luminosity(SED* sed) const
{
    // get the normalized mean specific luminosity for the SED convolved with the band
    Array lambdav, pv;
    sed->specificLuminosityArray(lambdav, pv, _band->wavelengthRange());
    double LlambdaSED = _band->meanSpecificLuminosity(lambdav, pv);
    if (LlambdaSED <= 0) throw FATALERROR("The normalization band is outside of the SED's wavelength range");

    // convert the user-configured specific luminosity to per-wavelength units
    double lambda = _band->pivotWavelength();
    double LlambdaUser = 0.;
    switch (_unitStyle)
    {
    case UnitStyle::wavelengthmonluminosity:
        LlambdaUser = _specificLuminosity;
        break;
    case UnitStyle::frequencymonluminosity:
        LlambdaUser = _specificLuminosity * Constants::c() / lambda / lambda;
        break;
    case UnitStyle::neutralmonluminosity:
        LlambdaUser = _specificLuminosity / lambda;
        break;
    }

    // return the ratio
    return LlambdaUser / LlambdaSED;
}

//////////////////////////////////////////////////////////////////////
