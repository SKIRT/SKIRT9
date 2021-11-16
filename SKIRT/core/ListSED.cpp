/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "ListSED.hpp"
#include "FatalError.hpp"
#include "NR.hpp"

//////////////////////////////////////////////////////////////////////

void ListSED::getWavelengthsAndLuminosities(Array& lambdav, Array& pv) const
{
    // verify number of configured parameters
    if (_wavelengths.size() != _specificLuminosities.size())
        throw FATALERROR("Number of listed luminosities does not match number of listed wavelengths");

    // copy the wavelengths and specific luminosities from the configuration
    NR::assign(lambdav, _wavelengths);
    NR::assign(pv, _specificLuminosities);

    // convert the user-configured specific luminosity to per-wavelength units without worrying about constant scale
    switch (_unitStyle)
    {
        case UnitStyle::neutralmonluminosity: pv /= lambdav; break;
        case UnitStyle::wavelengthmonluminosity: break;
        case UnitStyle::frequencymonluminosity: pv /= (lambdav * lambdav); break;
        case UnitStyle::energymonluminosity: pv /= (lambdav * lambdav * lambdav); break;
    }
}

//////////////////////////////////////////////////////////////////////
