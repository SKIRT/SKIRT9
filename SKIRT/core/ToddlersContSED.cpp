/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "ToddlersContSED.hpp"
#include "ToddlersContSEDFamily.hpp"
#include "NR.hpp"

////////////////////////////////////////////////////////////////////

const SEDFamily* ToddlersContSED::getFamilyAndParameters(Array& parameters)
{
    // set the parameters using arbitrary scaling
    NR::assign(parameters, _age, _metallicity, _SFE, _cloudNumDensity, 1.);

    // construct the library of SED models
    return new ToddlersContSEDFamily(this);
}

////////////////////////////////////////////////////////////////////
