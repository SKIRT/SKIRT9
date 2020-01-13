/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "BruzualCharlotSED.hpp"
#include "BruzualCharlotSEDFamily.hpp"
#include "NR.hpp"

////////////////////////////////////////////////////////////////////

const SEDFamily* BruzualCharlotSED::getFamilyAndParameters(Array& parameters)
{
    // set the parameters using arbitrary scaling
    NR::assign(parameters, 1., _metallicity, _age);

    // construct the library of SED models
    return new BruzualCharlotSEDFamily(
        this, _imf == IMF::Chabrier ? BruzualCharlotSEDFamily::IMF::Chabrier : BruzualCharlotSEDFamily::IMF::Salpeter,
        BruzualCharlotSEDFamily::Resolution::High);
}

////////////////////////////////////////////////////////////////////
