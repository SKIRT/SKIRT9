/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "MagneticFieldProbe.hpp"
#include "Configuration.hpp"
#include "MediumSystem.hpp"
#include "ProbeFormBridge.hpp"

////////////////////////////////////////////////////////////////////

void MagneticFieldProbe::probe()
{
    if (find<Configuration>()->hasMagneticField())
    {
        // locate the medium system
        auto ms = find<MediumSystem>();

        // define the call-back functions that will be passed to the bridge
        auto valueInCell = [ms](int m) { return ms->magneticField(m); };
        auto weightInCell = [ms](int m) { return ms->massDensity(m); };

        // construct a bridge and tell it to produce output
        ProbeFormBridge bridge(this, form());
        bridge.writeQuantity("B", "magneticfield", "magnetic field", "density-weighted magnetic field", valueInCell,
                             weightInCell);
    }
}

////////////////////////////////////////////////////////////////////
