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
        // define a call-back function to retrieve the density for a given entity
        auto densityInEntity = [snapshot, this](int m) {
            double mass = massType() == MassType::CurrentMass ? snapshot->currentMass(m) : snapshot->initialMass(m);
            return mass / snapshot->volume(m);
        };

        // construct a bridge and produce output
        ProbeFormBridge bridge(this, form());
        string stype = massType() == MassType::CurrentMass ? "current" : "initial";
        bridge.writeQuantity(sh + "_rho", sh + "_Sigma", "massvolumedensity", "masssurfacedensity",
                             stype + " mass density", stype + " mass surface density", snapshot, densityInEntity);
    }
}

////////////////////////////////////////////////////////////////////
