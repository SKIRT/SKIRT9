/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "ToddlersLineSED.hpp"
#include "ToddlersLineSEDFamily.hpp"
#include "NR.hpp"

////////////////////////////////////////////////////////////////////

const SEDFamily* ToddlersLineSED::getFamilyAndParameters(Array& parameters)
{
    // set the parameters using arbitrary scaling
    NR::assign(parameters, _age, _metallicity, _SFE, _cloudNumDensity, 1.);

    // construct the library of SED models
    return new ToddlersLineSEDFamily(this);
}

////////////////////////////////////////////////////////////////////
