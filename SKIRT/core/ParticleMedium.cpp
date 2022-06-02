/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "ParticleMedium.hpp"
#include "ParticleSnapshot.hpp"

////////////////////////////////////////////////////////////////////

Snapshot* ParticleMedium::createAndOpenSnapshot()
{
    // create and open the snapshot
    auto snapshot = new ParticleSnapshot;
    snapshot->open(this, filename(), "smoothed particles");

    // honor custom column reordering
    snapshot->useColumns(useColumns());

    // configure the position and size columns
    snapshot->importPosition();
    snapshot->importSize();

    // configure the mass or number column
    switch (massType())
    {
        case MassType::Mass: snapshot->importMass(); break;
        case MassType::Number: snapshot->importNumber(); break;
    }

    // set the smoothing kernel
    snapshot->setSmoothingKernel(smoothingKernel());
    return snapshot;
}

////////////////////////////////////////////////////////////////////
