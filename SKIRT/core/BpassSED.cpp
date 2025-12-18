/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "BpassSED.hpp"
#include "BpassSEDFamily.hpp"
#include "NR.hpp"

////////////////////////////////////////////////////////////////////

const SEDFamily* BpassSED::getFamilyAndParameters(Array& parameters)
{
    // set the parameters using arbitrary scaling
    NR::assign(parameters, 1., _metallicity, _age);

    // construct the appropriate library of SED models
    auto imf = BpassSEDFamily::IMF::Chabrier100;
    switch (_imf)
    {
        case IMF::Chabrier100: imf = BpassSEDFamily::IMF::Chabrier100; break;
        case IMF::Chabrier300: imf = BpassSEDFamily::IMF::Chabrier300; break;
        case IMF::Kroupa100: imf = BpassSEDFamily::IMF::Kroupa100; break;
        case IMF::Kroupa300: imf = BpassSEDFamily::IMF::Kroupa300; break;
    }
    auto resolution = _resolution == Resolution::Downsampled ? BpassSEDFamily::Resolution::Downsampled
                                                             : BpassSEDFamily::Resolution::Original;

    return new BpassSEDFamily(this, imf, resolution);
}

////////////////////////////////////////////////////////////////////
