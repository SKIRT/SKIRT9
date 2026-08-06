/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "CellSource.hpp"
#include "CellSnapshot.hpp"

////////////////////////////////////////////////////////////////////

Snapshot* CellSource::createAndOpenSnapshot()
{
    // create and open the snapshot
    auto snapshot = new CellSnapshot;
    snapshot->open(this, filename(), "cuboidal source cells");

    // honor custom column reordering
    snapshot->useColumns(useColumns());

    // configure the cell bounding box columns
    snapshot->importBox();
    return snapshot;
}

////////////////////////////////////////////////////////////////////
