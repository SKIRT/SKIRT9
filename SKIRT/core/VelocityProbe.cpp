/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "VelocityProbe.hpp"
#include "Configuration.hpp"
#include "MediumSystem.hpp"
#include "ProbeFormBridge.hpp"
#include "TextOutFile.hpp"
#include "Units.hpp"

////////////////////////////////////////////////////////////////////

void VelocityProbe::probe()
{
    if (find<Configuration>()->hasMovingMedia())
    {
        // locate the medium system
        auto ms = find<MediumSystem>();

        // define the call-back functions that will be passed to the bridge
        auto addColumnDefinitions = [this](TextOutFile& outfile) {
            auto units = find<Units>();
            outfile.addColumn("x component of velocity", units->uvelocity());
            outfile.addColumn("y component of velocity", units->uvelocity());
            outfile.addColumn("z component of velocity", units->uvelocity());
        };
        auto valueInCell = [ms](int m) { return ms->bulkVelocity(m); };
        auto weightInCell = [ms](int m) { return ms->massDensity(m); };

        // construct a bridge and tell it to produce output
        ProbeFormBridge bridge(this, form());
        bridge.writeQuantity("v", "velocity", "medium velocity", "mass-weighted medium velocity", addColumnDefinitions,
                             valueInCell, weightInCell);
    }
}

////////////////////////////////////////////////////////////////////
