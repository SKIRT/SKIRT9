/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "ParticleGeometry.hpp"
#include "ParticleSnapshot.hpp"

////////////////////////////////////////////////////////////////////

Snapshot* ParticleGeometry::createAndOpenSnapshot()
{
    // create and open the snapshot
    auto snapshot = new ParticleSnapshot;
    snapshot->open(this, filename(), "smoothed particles");

    // honor custom column reordering
    snapshot->useColumns(useColumns());

    // configure the position, size and mass columns
    snapshot->importPosition();
    snapshot->importSize();
    snapshot->importMass();

    // set the smoothing kernel
    snapshot->setSmoothingKernel(smoothingKernel());
    return snapshot;
}

////////////////////////////////////////////////////////////////////
