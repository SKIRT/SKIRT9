/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "TetraMeshSource.hpp"
#include "TetraMeshSnapshot.hpp"

////////////////////////////////////////////////////////////////////

Snapshot* TetraMeshSource::createAndOpenSnapshot()
{
    // create and open the snapshot
    auto snapshot = new TetraMeshSnapshot;
    snapshot->open(this, filename(), "Tetra source sites");

    // honor custom column reordering
    snapshot->useColumns(useColumns());

    // configure the position columns
    snapshot->importPosition();

    // set the domain extent
    snapshot->setExtent(domain());
    return snapshot;
}

////////////////////////////////////////////////////////////////////
