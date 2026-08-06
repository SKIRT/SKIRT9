/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "VoronoiMeshGeometry.hpp"
#include "VoronoiMeshSnapshot.hpp"

////////////////////////////////////////////////////////////////////

Snapshot* VoronoiMeshGeometry::createAndOpenSnapshot()
{
    // create and open the snapshot
    _voronoiMeshSnapshot = new VoronoiMeshSnapshot;
    _voronoiMeshSnapshot->open(this, filename(), "Voronoi sites");

    // honor custom column reordering
    _voronoiMeshSnapshot->useColumns(useColumns());

    // configure the position columns
    _voronoiMeshSnapshot->importPosition();

    // configure the mass or density column
    switch (massType())
    {
        case MassType::MassDensity: _voronoiMeshSnapshot->importMassDensity(); break;
        case MassType::Mass: _voronoiMeshSnapshot->importMass(); break;
        case MassType::NumberDensity: _voronoiMeshSnapshot->importNumberDensity(); break;
        case MassType::Number: _voronoiMeshSnapshot->importNumber(); break;
    }

    // set the domain extent
    _voronoiMeshSnapshot->setExtent(domain());
    return _voronoiMeshSnapshot;
}

////////////////////////////////////////////////////////////////////

VoronoiMeshSnapshot* VoronoiMeshGeometry::voronoiMesh() const
{
    return _voronoiMeshSnapshot;
}

////////////////////////////////////////////////////////////////////
