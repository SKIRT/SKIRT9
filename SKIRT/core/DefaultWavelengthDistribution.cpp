/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "DefaultWavelengthDistribution.hpp"
#include "FatalError.hpp"
#include "Random.hpp"
#include "SourceWavelengthRangeInterface.hpp"

//////////////////////////////////////////////////////////////////////

void DefaultWavelengthDistribution::setupSelfBefore()
{
    WavelengthDistribution::setupSelfBefore();

    _range = interface<SourceWavelengthRangeInterface>()->wavelengthRange();
    if (_range.empty()) throw FATALERROR("Source wavelength range is empty");

    _logMin = log(_range.min());
    _logWidth = log(_range.max()) - log(_range.min());
}

//////////////////////////////////////////////////////////////////////

double DefaultWavelengthDistribution::probability(double wavelength) const
{
    if (_range.containsFuzzy(wavelength))
        return 1. / (_logWidth * wavelength);
    else
        return 0.;
}

//////////////////////////////////////////////////////////////////////

double DefaultWavelengthDistribution::generateWavelength() const
{
    return exp(_logMin + _logWidth * random()->uniform());
}

//////////////////////////////////////////////////////////////////////
