/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "AdaptiveMeshSource.hpp"
#include "AdaptiveMeshSnapshot.hpp"

////////////////////////////////////////////////////////////////////

Snapshot* AdaptiveMeshSource::createAndOpenSnapshot()
{
    // create and open the snapshot
    auto snapshot = new AdaptiveMeshSnapshot;
    snapshot->open(this, filename(), "adaptive mesh source cells");

    // honor custom column reordering
    snapshot->useColumns(useColumns());

    // set the domain extent
    snapshot->setExtent(domain());
    return snapshot;
}

////////////////////////////////////////////////////////////////////
