/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "ImportedMediumDensityProbe.hpp"
#include "EntityCollection.hpp"
#include "ProbeFormBridge.hpp"
#include "Snapshot.hpp"

////////////////////////////////////////////////////////////////////

void ImportedMediumDensityProbe::probeImportedMedium(string sh, const ImportedMedium* /*medium*/,
                                                     const Snapshot* snapshot)
{
    // construct a bridge and produce output
    ProbeFormBridge bridge(this, form());
    if (snapshot->holdsNumber())
    {
        bridge.writeQuantity(sh + "_n", sh + "_N", "numbervolumedensity", "numbersurfacedensity", "number density",
                             "column density", vector<const Snapshot*>{snapshot},
                             [](const Snapshot* s, int m) { return s->density(m); });
    }
    else
    {
        bridge.writeQuantity(sh + "_rho", sh + "_Sigma", "massvolumedensity", "masssurfacedensity", "mass density",
                             "mass surface density", vector<const Snapshot*>{snapshot},
                             [](const Snapshot* s, int m) { return s->density(m); });
    }
}

////////////////////////////////////////////////////////////////////
