/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "AdaptiveMeshGeometry.hpp"
#include "AdaptiveMeshSnapshot.hpp"

////////////////////////////////////////////////////////////////////

Snapshot* AdaptiveMeshGeometry::createAndOpenSnapshot()
{
    // create and open the snapshot
    auto snapshot = new AdaptiveMeshSnapshot;
    snapshot->open(this, filename(), "AMR cells");

    // configure the mass or density column
    if (useMass()) snapshot->importMass();
    else snapshot->importDensity();

    // set the domain extent
    snapshot->setExtent(domain());
    return snapshot;
}

////////////////////////////////////////////////////////////////////
