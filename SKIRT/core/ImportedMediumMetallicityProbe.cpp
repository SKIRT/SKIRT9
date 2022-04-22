/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "ImportedMediumMetallicityProbe.hpp"
#include "EntityCollection.hpp"
#include "ImportedMedium.hpp"
#include "ProbeFormBridge.hpp"
#include "Snapshot.hpp"

////////////////////////////////////////////////////////////////////

void ImportedMediumMetallicityProbe::probeImportedMedium(string sh, const ImportedMedium* medium,
                                                         const Snapshot* snapshot)
{
    if (snapshot->hasMetallicity())
    {
        // construct a bridge
        ProbeFormBridge bridge(this, form());

        // define call-back functions to retrieve the probed quantity and corresponding weight for a given entity
        // for dust media, use the gas density rather than the dust density for weighting the probed quantity
        bool dust = medium->mix()->isDust();
        auto getValue = [snapshot](int m) { return snapshot->metallicity(m); };
        auto getWeight = [snapshot, dust](int m) {
            double w = snapshot->density(m);
            if (dust) w /= snapshot->metallicity(m);
            return w;
        };

        // define the call-back function to retrieve an averaged quantity value at a given position
        auto valueAtPosition = [snapshot, getValue, getWeight](Position bfr) {
            thread_local EntityCollection entities;  // can be reused for all queries in a given execution thread
            snapshot->getEntities(entities, bfr);
            return entities.average(getValue, getWeight);
        };

        // define the call-back function to retrieve an averaged quantity value along a given path
        auto valueAlongPath = [snapshot, getValue, getWeight](Position bfr, Direction bfk) {
            thread_local EntityCollection entities;  // can be reused for all queries in a given execution thread
            snapshot->getEntities(entities, bfr, bfk);
            return entities.average(getValue, getWeight);
        };

        // produce output
        bridge.writeQuantity(sh + "_n", sh + "_N", "dimensionless", "dimensionless", "metallicity",
                             "density-weighted metallicity", valueAtPosition, valueAlongPath);
    }
}

////////////////////////////////////////////////////////////////////
