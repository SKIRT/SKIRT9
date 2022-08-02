/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "ImportedSourceMetallicityProbe.hpp"
#include "EntityCollection.hpp"
#include "ProbeFormBridge.hpp"
#include "Snapshot.hpp"

////////////////////////////////////////////////////////////////////

void ImportedSourceMetallicityProbe::probeImportedSourceWeighted(
    string sweight, const vector<const Snapshot*>& snapshots,
    std::function<double(const Snapshot* snapshot, int m)> weight)
{
    // verify that all snapshots offer the metallicity property
    bool haveMetallicity = true;
    for (auto snapshot : snapshots) haveMetallicity &= snapshot->hasMetallicity();

    if (haveMetallicity)
    {
        // construct a bridge and produce output
        ProbeFormBridge bridge(this, form());
        bridge.writeQuantity(
            "Z", "dimensionless", "metallicity", sweight + "-weighted metallicity", snapshots,
            [](const Snapshot* snapshot, int m) { return snapshot->metallicity(m); }, weight);
    }
}

////////////////////////////////////////////////////////////////////
