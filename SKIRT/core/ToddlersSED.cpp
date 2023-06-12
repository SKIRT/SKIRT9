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
    NR::assign(parameters, _age, _metallicity, _SFE, _cloudNumDensity, 1.);

    // construct the library of SED models
    return new ToddlersSEDFamily(
        this,
        _pahfraction == PAHFraction::High ? ToddlersSEDFamily::PAHFraction::High : ToddlersSEDFamily::PAHFraction::Low,
        _resolution == Resolution::Low ? ToddlersSEDFamily::Resolution::Low : ToddlersSEDFamily::Resolution::High);
}

////////////////////////////////////////////////////////////////////
