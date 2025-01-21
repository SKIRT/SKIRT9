/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "CylindricalCellSource.hpp"
#include "CylindricalCellSnapshot.hpp"

////////////////////////////////////////////////////////////////////

Snapshot* CylindricalCellSource::createAndOpenSnapshot()
{
    // create and open the snapshot
    auto snapshot = new CylindricalCellSnapshot;
    snapshot->open(this, filename(), "cylindrical source cells");

    // honor custom column reordering
    snapshot->useColumns(useColumns());

    // configure the cell bounding box columns
    snapshot->setCoordinateSystem(CylindricalCellSnapshot::CoordinateSystem::CYLINDRICAL);
    snapshot->importBox();
    if (autoRevolve()) snapshot->setNumAutoRevolveBins(numAutoRevolveBins());

    return snapshot;
}

////////////////////////////////////////////////////////////////////
