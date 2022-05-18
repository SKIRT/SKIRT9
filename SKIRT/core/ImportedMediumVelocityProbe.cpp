/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "ImportedMediumVelocityProbe.hpp"
#include "EntityCollection.hpp"
#include "ImportedMedium.hpp"
#include "ProbeFormBridge.hpp"
#include "Snapshot.hpp"

////////////////////////////////////////////////////////////////////

void ImportedMediumVelocityProbe::probeImportedMedium(string sh, const ImportedMedium* medium, const Snapshot* snapshot)
{
    if (snapshot->hasVelocity())
    {
        // define call-back functions to retrieve the probed quantity and corresponding weight for a given entity
        // for dust media, use the gas density rather than the dust density for weighting the probed quantity
        bool dust = medium->mix()->isDust() && snapshot->hasMetallicity();
        auto getValue = [snapshot](int m) { return snapshot->velocity(m); };
        auto getWeight = [snapshot, dust](int m) {
            double w = snapshot->density(m);
            if (dust) w /= snapshot->metallicity(m);
            return w;
        };

        // construct a bridge and produce output
        ProbeFormBridge bridge(this, form());
        bridge.writeQuantity(sh + "_v", "velocity", "velocity", "density-weighted velocity", snapshot, getValue,
                             getWeight);
    }
}

////////////////////////////////////////////////////////////////////
