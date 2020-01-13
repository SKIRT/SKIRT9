/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "ListWavelengthDistribution.hpp"
#include "FatalError.hpp"
#include "NR.hpp"

//////////////////////////////////////////////////////////////////////

void ListWavelengthDistribution::getWavelengthsAndProbabilities(Array& lambdav, Array& pv) const
{
    // verify number of configured parameters
    if (_wavelengths.size() != _probabilities.size())
        throw FATALERROR("Number of listed probabilities does not match number of listed wavelengths");

    // copy the wavelengths and probabilities from the configuration
    NR::assign(lambdav, _wavelengths);
    NR::assign(pv, _probabilities);

    // convert the user-configured probabilities to per-wavelength units without worrying about constant scale
    switch (_unitStyle)
    {
        case UnitStyle::wavelengthmonluminosity: break;
        case UnitStyle::frequencymonluminosity: pv /= (lambdav * lambdav); break;
        case UnitStyle::neutralmonluminosity: pv /= lambdav; break;
    }
}

//////////////////////////////////////////////////////////////////////
