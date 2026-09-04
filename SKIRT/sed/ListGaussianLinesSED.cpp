/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "ListGaussianLinesSED.hpp"
#include "NR.hpp"

//////////////////////////////////////////////////////////////////////

void ListGaussianLinesSED::getLineProperties(Array& lambdav, Array& dispersionv, Array& Lv) const
{
    // copy the wavelengths, dispersions and luminosities from the configuration
    NR::assign(lambdav, _wavelengths);
    NR::assign(dispersionv, _dispersions);
    NR::assign(Lv, _luminosities);
}

//////////////////////////////////////////////////////////////////////
