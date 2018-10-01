/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "RangeWavelengthDistribution.hpp"
#include "Configuration.hpp"
#include "FatalError.hpp"

//////////////////////////////////////////////////////////////////////

void RangeWavelengthDistribution::setupSelfBefore()
{
    WavelengthDistribution::setupSelfBefore();

    _range.set(_minWavelength, _maxWavelength);
    _range.intersect(find<Configuration>()->sourceWavelengthRange());
    if (_range.empty()) throw FATALERROR("Wavelength distribution range does not overlap source wavelength range");
}

//////////////////////////////////////////////////////////////////////
