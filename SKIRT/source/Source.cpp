/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "Source.hpp"
#include "Configuration.hpp"
#include "Log.hpp"
#include "Random.hpp"
#include "StringUtils.hpp"
#include "Units.hpp"

//////////////////////////////////////////////////////////////////////

void Source::setupSelfBefore()
{
    SimulationItem::setupSelfBefore();

    _random = find<Random>();
}

//////////////////////////////////////////////////////////////////////

void Source::informAvailableWavelengthRange(Range available, string itemType)
{
    auto configured = find<Configuration>()->sourceWavelengthRange();
    const double fuzzy = 0.01;  // 1%
    if (!available.containsFuzzy(configured.min(), fuzzy) || !available.containsFuzzy(configured.max(), fuzzy))
    {
        auto units = find<Units>();
        auto outstring = [units](double wavelength) {
            return StringUtils::toString(units->owavelength(wavelength), 'g', 3) + " " + units->uwavelength();
        };

        auto log = find<Log>();
        log->warning(itemType + " available spectrum does not fully cover the configured source wavelength range");
        log->warning("  Configured: " + outstring(configured.min()) + " - " + outstring(configured.max()));
        log->warning("  Available: " + outstring(available.min()) + " - " + outstring(available.max()));
    }
}

//////////////////////////////////////////////////////////////////////

void Source::prepareForLaunch(double /*sourceBias*/, size_t /*firstIndex*/, size_t /*numIndices*/) {}

//////////////////////////////////////////////////////////////////////
