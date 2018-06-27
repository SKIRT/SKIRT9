/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "ParticleSource.hpp"
#include "ParticleSnapshot.hpp"

////////////////////////////////////////////////////////////////////

Snapshot* ParticleSource::createAndOpenSnapshot()
{
    // create and open the snapshot
    auto snapshot = new ParticleSnapshot;
    snapshot->open(this, filename(), "smoothed source particles");

    // set the smoothing kernel
    snapshot->setSmoothingKernel(smoothingKernel());
    return snapshot;
}

////////////////////////////////////////////////////////////////////
