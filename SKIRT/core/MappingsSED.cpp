/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "MappingsSED.hpp"
#include "MappingsSEDFamily.hpp"
#include "NR.hpp"

////////////////////////////////////////////////////////////////////

const SEDFamily* MappingsSED::getFamilyAndParameters(Array& parameters)
{
    // set the parameters using arbitrary scaling
    NR::assign(parameters, 1., _metallicity, _compactness, _pressure, _coveringFactor);

    // construct the library of SED models
    return new MappingsSEDFamily(this);
}

////////////////////////////////////////////////////////////////////
