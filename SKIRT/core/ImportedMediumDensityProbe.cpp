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
    // construct a bridge
    ProbeFormBridge bridge(this, form());

    // define the call-back function to retrieve a density value at a given position
    auto valueAtPosition = [snapshot](Position bfr) {
        thread_local EntityCollection entities;  // can be reused for all queries in a given execution thread
        snapshot->getEntities(entities, bfr);
        return entities.accumulate([snapshot](int m) { return snapshot->density(m); });
    };

    // define the call-back function to retrieve a surface density value along a given path
    auto valueAlongPath = [snapshot](Position bfr, Direction bfk) {
        thread_local EntityCollection entities;  // can be reused for all queries in a given execution thread
        snapshot->getEntities(entities, bfr, bfk);
        return entities.accumulate([snapshot](int m) { return snapshot->density(m); });
    };

    // produce output
    if (snapshot->holdsNumber())
    {
        bridge.writeQuantity(sh + "_n", sh + "_N", "numbervolumedensity", "numbersurfacedensity", "number density",
                             "column density", valueAtPosition, valueAlongPath);
    }
    else
    {
        bridge.writeQuantity(sh + "_rho", sh + "_Sigma", "massvolumedensity", "masssurfacedensity", "mass density",
                             "mass surface density", valueAtPosition, valueAlongPath);
    }
}

////////////////////////////////////////////////////////////////////
