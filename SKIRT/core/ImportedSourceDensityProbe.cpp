/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "ImportedSourceDensityProbe.hpp"
#include "EntityCollection.hpp"
#include "ProbeFormBridge.hpp"
#include "Snapshot.hpp"

////////////////////////////////////////////////////////////////////

void ImportedSourceDensityProbe::probeImportedSource(string sh, const ImportedSource* /*source*/,
                                                     const Snapshot* snapshot)
{
    if ((massType() == MassType::InitialMass && snapshot->hasInitialMass())
        || (massType() == MassType::CurrentMass && snapshot->hasCurrentMass()))
    {
        // construct a bridge
        ProbeFormBridge bridge(this, form());

        // define a call-back function to retrieve the density for a given entity
        auto density = [snapshot, this](int m) {
            double mass = massType() == MassType::CurrentMass ? snapshot->currentMass(m) : snapshot->initialMass(m);
            return mass / snapshot->volume(m);
        };

        // define the call-back function to retrieve a density value at a given position
        auto valueAtPosition = [snapshot, density](Position bfr) {
            thread_local EntityCollection entities;  // can be reused for all queries in a given execution thread
            snapshot->getEntities(entities, bfr);
            return entities.accumulate(density);
        };

        // define the call-back function to retrieve a surface density value along a given path
        auto valueAlongPath = [snapshot, density](Position bfr, Direction bfk) {
            thread_local EntityCollection entities;  // can be reused for all queries in a given execution thread
            snapshot->getEntities(entities, bfr, bfk);
            return entities.accumulate(density);
        };

        // produce output
        string stype = massType() == MassType::CurrentMass ? "current" : "initial";
        bridge.writeQuantity(sh + "_rho", sh + "_Sigma", "massvolumedensity", "masssurfacedensity",
                             stype + " mass density", stype + " mass surface density", valueAtPosition, valueAlongPath);
    }
}

////////////////////////////////////////////////////////////////////
