/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "SphericalCellSource.hpp"
#include "SphericalCellSnapshot.hpp"

////////////////////////////////////////////////////////////////////

Snapshot* SphericalCellSource::createAndOpenSnapshot()
{
    // create and open the snapshot
    auto snapshot = new SphericalCellSnapshot;
    snapshot->open(this, filename(), "spherical source cells");

    // honor custom column reordering
    snapshot->useColumns(useColumns());

    // configure the cell bounding box columns
    snapshot->setCoordinateSystem(SphericalCellSnapshot::CoordinateSystem::SPHERICAL);
    snapshot->importBox();
    switch (autoRevolve())
    {
        case AutoRevolveType::None: break;
        case AutoRevolveType::Inclination: snapshot->setNumAutoRevolveBins(numInclinationRevolveBins(), 0); break;
        case AutoRevolveType::Azimuth: snapshot->setNumAutoRevolveBins(0, numAzimuthRevolveBins()); break;
        case AutoRevolveType::InclinationAndAzimuth:
            snapshot->setNumAutoRevolveBins(numInclinationRevolveBins(), numAzimuthRevolveBins());
            break;
    }

    return snapshot;
}

////////////////////////////////////////////////////////////////////
