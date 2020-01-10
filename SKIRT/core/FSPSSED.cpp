/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "FSPSSED.hpp"
#include "FSPSSEDFamily.hpp"
#include "NR.hpp"

////////////////////////////////////////////////////////////////////

const SEDFamily* FSPSSED::getFamilyAndParameters(Array& parameters)
{
    // set the parameters using arbitrary scaling
    NR::assign(parameters, 1., _metallicity, _age);

    // translate IMF enum
    FSPSSEDFamily::IMF imf = FSPSSEDFamily::IMF::Chabrier;
    switch (_imf)
    {
        case IMF::Chabrier: imf = FSPSSEDFamily::IMF::Chabrier; break;
        case IMF::Kroupa: imf = FSPSSEDFamily::IMF::Kroupa; break;
        case IMF::Salpeter: imf = FSPSSEDFamily::IMF::Salpeter; break;
    }

    // construct the library of SED models
    return new FSPSSEDFamily(this, imf);
}

////////////////////////////////////////////////////////////////////
