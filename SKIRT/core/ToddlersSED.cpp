/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "ToddlersSED.hpp"
#include "NR.hpp"
#include "ToddlersSEDFamily.hpp"

////////////////////////////////////////////////////////////////////

const SEDFamily* ToddlersSED::getFamilyAndParameters(Array& parameters)
{
    // set the parameters using arbitrary scaling
    NR::assign(parameters, _metallicity, _SFE, _cloudNumDensity, 1.);

    // construct the library of SED models
    return new ToddlersSEDFamily(this, ToddlersSEDFamily::SedMode::SFRNormalized,
                                 ToddlersSEDFamily::StellarTemplate::SB99Kroupa100Sin, true,
                                 ToddlersSEDFamily::Resolution::High, ToddlersSEDFamily::SFRPeriod::Period10Myr);
}

////////////////////////////////////////////////////////////////////
