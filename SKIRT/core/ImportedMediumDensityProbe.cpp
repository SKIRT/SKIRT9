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
                             "column density", snapshot, [snapshot](int m) { return snapshot->density(m); });
    }
    else
    {
        bridge.writeQuantity(sh + "_rho", sh + "_Sigma", "massvolumedensity", "masssurfacedensity", "mass density",
                             "mass surface density", snapshot, [snapshot](int m) { return snapshot->density(m); });
    }
}

////////////////////////////////////////////////////////////////////
