/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "SED.hpp"
#include "Configuration.hpp"
#include "FatalError.hpp"
#include "Random.hpp"

////////////////////////////////////////////////////////////////////

void SED::setupSelfBefore()
{
    SimulationItem::setupSelfBefore();

    _random = find<Random>();
}

////////////////////////////////////////////////////////////////////

Range SED::normalizationWavelengthRange() const
{
    Range range = find<Configuration>()->sourceWavelengthRange();
    range.intersect(intrinsicWavelengthRange());
    if (range.empty()) throw FATALERROR("Intrinsic SED wavelength range does not overlap source wavelength range");
    return range;
}

////////////////////////////////////////////////////////////////////
