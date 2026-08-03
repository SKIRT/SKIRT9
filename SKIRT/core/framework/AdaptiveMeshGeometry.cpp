/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "AdaptiveMeshGeometry.hpp"
#include "AdaptiveMeshSnapshot.hpp"

////////////////////////////////////////////////////////////////////

Snapshot* AdaptiveMeshGeometry::createAndOpenSnapshot()
{
    // create and open the snapshot
    _adaptiveMeshSnapshot = new AdaptiveMeshSnapshot;
    _adaptiveMeshSnapshot->open(this, filename(), "adaptive mesh cells");

    // honor custom column reordering
    _adaptiveMeshSnapshot->useColumns(useColumns());

    // configure the mass or density column
    switch (massType())
    {
        case MassType::MassDensity: _adaptiveMeshSnapshot->importMassDensity(); break;
        case MassType::Mass: _adaptiveMeshSnapshot->importMass(); break;
        case MassType::NumberDensity: _adaptiveMeshSnapshot->importNumberDensity(); break;
        case MassType::Number: _adaptiveMeshSnapshot->importNumber(); break;
    }

    // set the domain extent
    _adaptiveMeshSnapshot->setExtent(domain());
    return _adaptiveMeshSnapshot;
}

////////////////////////////////////////////////////////////////////

AdaptiveMeshSnapshot* AdaptiveMeshGeometry::adaptiveMesh() const
{
    return _adaptiveMeshSnapshot;
}

////////////////////////////////////////////////////////////////////
