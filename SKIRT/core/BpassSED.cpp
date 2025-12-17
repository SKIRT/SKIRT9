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

    // construct the library of SED models
    if (_resolution == Resolution::Downsampled)
    {
        switch (_imf)
        {
            case IMF::Chabrier100:
                return new BpassSEDFamily(this,BpassSEDFamily::IMF::Chabrier100,
                                          BpassSEDFamily::Resolution::Downsampled);
            case IMF::Chabrier300:
                return new BpassSEDFamily(this,BpassSEDFamily::IMF::Chabrier300,
                                          BpassSEDFamily::Resolution::Downsampled);
            case IMF::Kroupa100:
                return new BpassSEDFamily(this,BpassSEDFamily::IMF::Kroupa100,
                                          BpassSEDFamily::Resolution::Downsampled);
            case IMF::Kroupa300:
                return new BpassSEDFamily(this,BpassSEDFamily::IMF::Kroupa300,
                                          BpassSEDFamily::Resolution::Downsampled);
        }
    }
    else
    {
        switch (_imf)
        {
            case IMF::Chabrier100:
                return new BpassSEDFamily(this,BpassSEDFamily::IMF::Chabrier100,
                                          BpassSEDFamily::Resolution::Original);
            case IMF::Chabrier300:
                return new BpassSEDFamily(this,BpassSEDFamily::IMF::Chabrier300,
                                          BpassSEDFamily::Resolution::Original);
            case IMF::Kroupa100:
                return new BpassSEDFamily(this,BpassSEDFamily::IMF::Kroupa100,
                                          BpassSEDFamily::Resolution::Original);
            case IMF::Kroupa300:
                return new BpassSEDFamily(this,BpassSEDFamily::IMF::Kroupa300,
                                          BpassSEDFamily::Resolution::Original);
        }
    }
}

////////////////////////////////////////////////////////////////////
