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

    // honor custom column reordering
    snapshot->useColumns(useColumns());

    // configure the position and size columns
    snapshot->importPosition();
    snapshot->importSize();

    // set the smoothing kernel
    snapshot->setSmoothingKernel(smoothingKernel());
    return snapshot;
}

////////////////////////////////////////////////////////////////////
