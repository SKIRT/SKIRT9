/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "ImportedSourceVelocityProbe.hpp"
#include "EntityCollection.hpp"
#include "ProbeFormBridge.hpp"
#include "Snapshot.hpp"

////////////////////////////////////////////////////////////////////

void ImportedSourceVelocityProbe::probeImportedSourceWeighted(
    string sweight, const vector<const Snapshot*>& snapshots,
    std::function<double(const Snapshot* snapshot, int m)> weight)
{
    // verify that all snapshots offer a velocity
    bool haveVelocity = true;
    for (auto snapshot : snapshots) haveVelocity &= snapshot->hasVelocity();

    if (haveVelocity)
    {
        // construct a bridge and produce output
        ProbeFormBridge bridge(this, form());
        bridge.writeQuantity(
            "v", "velocity", "velocity", sweight + "-weighted velocity", snapshots,
            [](const Snapshot* snapshot, int m) { return snapshot->velocity(m); }, weight);
    }
}

////////////////////////////////////////////////////////////////////
