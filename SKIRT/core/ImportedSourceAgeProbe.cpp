/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "ImportedSourceAgeProbe.hpp"
#include "EntityCollection.hpp"
#include "ProbeFormBridge.hpp"
#include "Snapshot.hpp"

////////////////////////////////////////////////////////////////////

void ImportedSourceAgeProbe::probeImportedSourceWeighted(string sh, string sweight, const Snapshot* snapshot,
                                                         std::function<double(int)> weight)
{
    if (snapshot->hasAge())
    {
        // construct a bridge and produce output
        ProbeFormBridge bridge(this, form());
        bridge.writeQuantity(
            sh + "_t", "time", "age", sweight + "-weighted age", snapshot,
            [snapshot](int m) { return snapshot->age(m); }, weight);
    }
}

////////////////////////////////////////////////////////////////////
