/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "ImportedSourceMetallicityProbe.hpp"
#include "EntityCollection.hpp"
#include "ProbeFormBridge.hpp"
#include "Snapshot.hpp"

////////////////////////////////////////////////////////////////////

void ImportedSourceMetallicityProbe::probeImportedSourceWeighted(string sh, string sweight, const Snapshot* snapshot,
                                                                 std::function<double(int)> weight)
{
    if (snapshot->hasMetallicity())
    {
        // construct a bridge and produce output
        ProbeFormBridge bridge(this, form());
        bridge.writeQuantity(
            sh + "_Z", "dimensionless", "metallicity", sweight + "-weighted metallicity", snapshot,
            [snapshot](int m) { return snapshot->metallicity(m); }, weight);
    }
}

////////////////////////////////////////////////////////////////////
