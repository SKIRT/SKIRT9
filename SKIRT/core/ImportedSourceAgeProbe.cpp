/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "ImportedSourceAgeProbe.hpp"
#include "EntityCollection.hpp"
#include "ProbeFormBridge.hpp"
#include "Snapshot.hpp"

////////////////////////////////////////////////////////////////////

void ImportedSourceAgeProbe::probeImportedSourceWeighted(string sweight, const vector<const Snapshot*>& snapshots,
                                                         std::function<double(const Snapshot* snapshot, int m)> weight)
{
    // verify that all snapshots offer the age property
    bool haveAge = true;
    for (auto snapshot : snapshots) haveAge &= snapshot->hasAge();

    if (haveAge)
    {
        // construct a bridge and produce output
        ProbeFormBridge bridge(this, form());
        bridge.writeQuantity(
            "t", "time", "age", sweight + "-weighted age", snapshots,
            [](const Snapshot* snapshot, int m) { return snapshot->age(m); }, weight);
    }
}

////////////////////////////////////////////////////////////////////
