/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "BandLuminosityNormalization.hpp"
#include "ContSED.hpp"
#include "FatalError.hpp"
#include "Units.hpp"

////////////////////////////////////////////////////////////////////

double BandLuminosityNormalization::luminosityForSED(SED* sed) const
{
    auto contsed = dynamic_cast<ContSED*>(sed);
    if (!contsed) throw FATALERROR("Cannot use band luminosity normalization for a line emission spectrum");

    // get the normalized mean specific luminosity for the SED convolved with the band
    Array lambdav, pv;
    contsed->specificLuminosityArray(lambdav, pv, _band->wavelengthRange());
    double LlambdaSED = _band->meanSpecificLuminosity(lambdav, pv);
    if (LlambdaSED <= 0) throw FATALERROR("The normalization band is outside of the SED's wavelength range");

    // convert the user-configured specific luminosity to per-wavelength units
    double wavelength = _band->pivotWavelength();
    double LlambdaUser = Units::fromFluxStyle(wavelength, _specificLuminosity, Units::fluxStyle(_unitStyle));

    // return the ratio
    return LlambdaUser / LlambdaSED;
}

//////////////////////////////////////////////////////////////////////
