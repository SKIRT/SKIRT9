/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "Starburst99SED.hpp"
#include "NR.hpp"
#include "Starburst99SEDFamily.hpp"

////////////////////////////////////////////////////////////////////

const SEDFamily* Starburst99SED::getFamilyAndParameters(Array& parameters)
{
    // set the parameters using arbitrary scaling
    NR::assign(parameters, 1., _metallicity, _age);

    // construct the library of SED models
    return new Starburst99SEDFamily(this);
}

////////////////////////////////////////////////////////////////////
