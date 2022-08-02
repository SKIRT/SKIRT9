/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "ImportedSourceDensityProbe.hpp"
#include "EntityCollection.hpp"
#include "ProbeFormBridge.hpp"
#include "Snapshot.hpp"

////////////////////////////////////////////////////////////////////

void ImportedSourceDensityProbe::probeImportedSources(const vector<const ImportedSource*>& /*sources*/,
                                                      const vector<const Snapshot*>& snapshots)
{
    // verify that all snapshots offer the requested mass type
    bool haveMassType = true;
    for (auto snapshot : snapshots)
    {
        haveMassType &= (massType() == MassType::InitialMass && snapshot->hasInitialMass())
                        || (massType() == MassType::CurrentMass && snapshot->hasCurrentMass());
    }

    if (haveMassType)
    {
        // define a call-back function to retrieve the density for a given entity
        auto densityInEntity = [this](const Snapshot* snapshot, int m) {
            double mass = massType() == MassType::CurrentMass ? snapshot->currentMass(m) : snapshot->initialMass(m);
            return mass / snapshot->volume(m);
        };

        // construct a bridge and produce output
        ProbeFormBridge bridge(this, form());
        string stype = massType() == MassType::CurrentMass ? "current" : "initial";
        bridge.writeQuantity("rho", "Sigma", "massvolumedensity", "masssurfacedensity", stype + " mass density",
                             stype + " mass surface density", snapshots, densityInEntity);
    }
}

////////////////////////////////////////////////////////////////////
