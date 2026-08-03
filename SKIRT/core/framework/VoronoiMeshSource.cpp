/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "VoronoiMeshSource.hpp"
#include "VoronoiMeshSnapshot.hpp"

////////////////////////////////////////////////////////////////////

Snapshot* VoronoiMeshSource::createAndOpenSnapshot()
{
    // create and open the snapshot
    auto snapshot = new VoronoiMeshSnapshot;
    snapshot->open(this, filename(), "Voronoi source sites");

    // honor custom column reordering
    snapshot->useColumns(useColumns());

    // configure the position columns
    snapshot->importPosition();

    // set the domain extent
    snapshot->setExtent(domain());
    return snapshot;
}

////////////////////////////////////////////////////////////////////
