/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "VelocityProbe.hpp"
#include "Configuration.hpp"
#include "MediumSystem.hpp"
#include "ProbeFormBridge.hpp"

////////////////////////////////////////////////////////////////////

void VelocityProbe::probe()
{
    if (find<Configuration>()->hasMovingMedia())
    {
        // locate the medium system
        auto ms = find<MediumSystem>();

        // define the call-back functions that will be passed to the bridge
        auto valueInCell = [ms](int m) { return ms->bulkVelocity(m); };
        auto weightInCell = [ms](int m) { return ms->massDensity(m); };

        // construct a bridge and tell it to produce output
        ProbeFormBridge bridge(this, form());
        bridge.writeQuantity("v", "velocity", "medium velocity", "density-weighted medium velocity", valueInCell,
                             weightInCell);
    }
}

////////////////////////////////////////////////////////////////////
