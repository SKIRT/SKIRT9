/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "VoronoiMeshGeometry.hpp"
#include "VoronoiMeshSnapshot.hpp"

////////////////////////////////////////////////////////////////////

Snapshot* VoronoiMeshGeometry::createAndOpenSnapshot()
{
    // create and open the snapshot
    auto snapshot = new VoronoiMeshSnapshot;
    snapshot->open(this, filename(), "Voronoi sites");

    // configure the mass or density column (position columns are configured by the snapshot itself)
    if (useMass()) snapshot->importMass();
    else snapshot->importDensity();

    // set the domain extent
    snapshot->setExtent(domain());
    return snapshot;
}

////////////////////////////////////////////////////////////////////
