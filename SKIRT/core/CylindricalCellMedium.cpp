/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "CylindricalCellMedium.hpp"
#include "CylindricalCellSnapshot.hpp"

////////////////////////////////////////////////////////////////////

Snapshot* CylindricalCellMedium::createAndOpenSnapshot()
{
    // create and open the snapshot
    auto snapshot = new CylindricalCellSnapshot;
    snapshot->open(this, filename(), "cylindrical cells");

    // honor custom column reordering
    snapshot->useColumns(useColumns());

    // configure the cell bounding box columns
    snapshot->setCoordinateSystem(CylindricalCellSnapshot::CoordinateSystem::CYLINDRICAL);
    snapshot->importBox();
    if (autoRevolve()) snapshot->setNumAutoRevolveBins(numAutoRevolveBins());

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
