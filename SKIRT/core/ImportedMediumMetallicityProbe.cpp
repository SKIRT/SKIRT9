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
        // define call-back functions to retrieve the probed quantity and corresponding weight for a given entity
        // for dust media, use the gas density rather than the dust density for weighting the probed quantity
        bool dust = medium->mix()->isDust();
        auto getValue = [](const Snapshot* snapshot, int m) { return snapshot->metallicity(m); };
        auto getWeight = [dust](const Snapshot* snapshot, int m) {
            double w = snapshot->density(m);
            if (dust) w /= snapshot->metallicity(m);
            return w;
        };

        // construct a bridge and produce output
        ProbeFormBridge bridge(this, form());
        bridge.writeQuantity(sh + "_Z", "dimensionless", "metallicity", "density-weighted metallicity",
                             vector<const Snapshot*>{snapshot}, getValue, getWeight);
    }
}

////////////////////////////////////////////////////////////////////
