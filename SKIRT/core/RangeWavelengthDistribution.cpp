/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "RangeWavelengthDistribution.hpp"
#include "FatalError.hpp"
#include "WavelengthRangeInterface.hpp"

//////////////////////////////////////////////////////////////////////

void RangeWavelengthDistribution::setupSelfBefore()
{
    WavelengthDistribution::setupSelfBefore();

    _range.set(_minWavelength, _maxWavelength);
    _range.intersect(interface<WavelengthRangeInterface>()->wavelengthRange());
    if (_range.empty()) throw FATALERROR("Wavelength distribution range does not overlap source wavelength range");
}

//////////////////////////////////////////////////////////////////////
