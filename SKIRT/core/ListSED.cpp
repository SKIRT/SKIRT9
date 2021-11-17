/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "ListSED.hpp"
#include "FatalError.hpp"
#include "NR.hpp"
#include "Units.hpp"

//////////////////////////////////////////////////////////////////////

void ListSED::getWavelengthsAndLuminosities(Array& lambdav, Array& pv) const
{
    // verify number of configured parameters
    if (_wavelengths.size() != _specificLuminosities.size())
        throw FATALERROR("Number of listed luminosities does not match number of listed wavelengths");

    // copy the user-configured wavelengths
    NR::assign(lambdav, _wavelengths);

    // convert the user-configured specific luminosities to per-wavelength units
    pv = Units::fromFluxStyle(lambdav, NR::array(_specificLuminosities), Units::fluxStyle(_unitStyle));
}

//////////////////////////////////////////////////////////////////////
