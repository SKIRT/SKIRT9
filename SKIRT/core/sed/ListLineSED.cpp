/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "ListLineSED.hpp"
#include "FatalError.hpp"
#include "NR.hpp"

//////////////////////////////////////////////////////////////////////

void ListLineSED::getWavelengthsAndLuminosities(Array& lambdav, Array& Lv) const
{
    // verify number of configured parameters
    if (_wavelengths.size() != _luminosities.size())
        throw FATALERROR("Number of listed luminosities does not match number of listed wavelengths");

    // copy the wavelengths and luminosities from the configuration
    NR::assign(lambdav, _wavelengths);
    NR::assign(Lv, _luminosities);
}

//////////////////////////////////////////////////////////////////////
