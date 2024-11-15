/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "TetraMeshGeometry.hpp"
#include "TetraMeshSnapshot.hpp"

////////////////////////////////////////////////////////////////////

Snapshot* TetraMeshGeometry::createAndOpenSnapshot()
{
    // create and open the snapshot
    _tetraMeshSnapshot = new TetraMeshSnapshot;
    _tetraMeshSnapshot->open(this, filename(), "Tetra sites");

    // honor custom column reordering
    _tetraMeshSnapshot->useColumns(useColumns());

    // configure the position columns
    _tetraMeshSnapshot->importPosition();

    // configure the mass or density column
    switch (massType())
    {
        case MassType::MassDensity: _tetraMeshSnapshot->importMassDensity(); break;
        case MassType::Mass: _tetraMeshSnapshot->importMass(); break;
        case MassType::NumberDensity: _tetraMeshSnapshot->importNumberDensity(); break;
        case MassType::Number: _tetraMeshSnapshot->importNumber(); break;
    }

    // set the domain extent
    _tetraMeshSnapshot->setExtent(domain());
    return _tetraMeshSnapshot;
}

////////////////////////////////////////////////////////////////////

TetraMeshSnapshot* TetraMeshGeometry::tetraMesh() const
{
    return _tetraMeshSnapshot;
}

////////////////////////////////////////////////////////////////////
