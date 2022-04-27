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
        // construct a bridge
        ProbeFormBridge bridge(this, form());

        // define the call-back function to retrieve an averaged quantity value at a given position
        auto valueAtPosition = [snapshot, weight](Position bfr) {
            thread_local EntityCollection entities;  // can be reused for all queries in a given execution thread
            snapshot->getEntities(entities, bfr);
            return entities.average([snapshot](int m) { return snapshot->velocity(m); }, weight);
        };

        // define the call-back function to retrieve an averaged quantity value along a given path
        auto valueAlongPath = [snapshot, weight](Position bfr, Direction bfk) {
            thread_local EntityCollection entities;  // can be reused for all queries in a given execution thread
            snapshot->getEntities(entities, bfr, bfk);
            return entities.average([snapshot](int m) { return snapshot->velocity(m); }, weight);
        };

        // produce output
        bridge.writeQuantity(sh + "_v", sh + "_v", "velocity", "velocity", "velocity", sweight + "-weighted velocity",
                             valueAtPosition, valueAlongPath);
    }
}

////////////////////////////////////////////////////////////////////
