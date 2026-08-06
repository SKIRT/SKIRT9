/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "ListWavelengthDistribution.hpp"
#include "FatalError.hpp"
#include "NR.hpp"
#include "Units.hpp"

//////////////////////////////////////////////////////////////////////

void ListWavelengthDistribution::getWavelengthsAndProbabilities(Array& lambdav, Array& pv) const
{
    // verify number of configured parameters
    if (_wavelengths.size() != _probabilities.size())
        throw FATALERROR("Number of listed probabilities does not match number of listed wavelengths");

    // copy the user-configured wavelengths
    NR::assign(lambdav, _wavelengths);

    // convert the user-configured probabilities to per-wavelength units
    pv = Units::fromFluxStyle(lambdav, NR::array(_probabilities), Units::fluxStyle(_unitStyle));
}

//////////////////////////////////////////////////////////////////////
