/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "ImportedSourceVelocityProbe.hpp"
#include "EntityCollection.hpp"
#include "ProbeFormBridge.hpp"
#include "Snapshot.hpp"

////////////////////////////////////////////////////////////////////

void ImportedSourceVelocityProbe::probeImportedSourceWeighted(string sh, string sweight, const Snapshot* snapshot,
                                                              std::function<double(int)> weight)
{
    if (snapshot->hasVelocity())
    {
        // construct a bridge and produce output
        ProbeFormBridge bridge(this, form());
        bridge.writeQuantity(
            sh + "_v", "velocity", "velocity", sweight + "-weighted velocity", snapshot,
            [snapshot](int m) { return snapshot->velocity(m); }, weight);
    }
}

////////////////////////////////////////////////////////////////////
