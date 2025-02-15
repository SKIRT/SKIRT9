/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "SphericalCellGeometry.hpp"
#include "SphericalCellSnapshot.hpp"

////////////////////////////////////////////////////////////////////

Snapshot* SphericalCellGeometry::createAndOpenSnapshot()
{
    // create and open the snapshot
    auto snapshot = new SphericalCellSnapshot;
    snapshot->open(this, filename(), "spherical cells");

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

    // configure the mass or density column
    switch (massType())
    {
        case MassType::MassDensity: snapshot->importMassDensity(); break;
        case MassType::Mass: snapshot->importMass(); break;
        case MassType::NumberDensity: snapshot->importNumberDensity(); break;
        case MassType::Number: snapshot->importNumber(); break;
    }
    return snapshot;
}

////////////////////////////////////////////////////////////////////
