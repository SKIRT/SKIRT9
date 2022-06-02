/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "VoronoiMeshMedium.hpp"
#include "Configuration.hpp"
#include "VoronoiMeshSnapshot.hpp"

////////////////////////////////////////////////////////////////////

Snapshot* VoronoiMeshMedium::createAndOpenSnapshot()
{
    // create and open the snapshot
    _voronoiMeshSnapshot = new VoronoiMeshSnapshot;
    _voronoiMeshSnapshot->open(this, filename(), "Voronoi sites");

    // honor custom column reordering
    _voronoiMeshSnapshot->useColumns(useColumns());

    // configure the position columns
    _voronoiMeshSnapshot->importPosition();

    // configure the density and/or mass column(s)
    bool bothDensityAndMass = false;
    switch (massType())
    {
        case MassType::MassDensity: _voronoiMeshSnapshot->importMassDensity(); break;
        case MassType::Mass: _voronoiMeshSnapshot->importMass(); break;
        case MassType::MassDensityAndMass:
            _voronoiMeshSnapshot->importMassDensity();
            _voronoiMeshSnapshot->importMass();
            bothDensityAndMass = true;
            break;

        case MassType::NumberDensity: _voronoiMeshSnapshot->importNumberDensity(); break;
        case MassType::Number: _voronoiMeshSnapshot->importNumber(); break;
        case MassType::NumberDensityAndNumber:
            _voronoiMeshSnapshot->importNumberDensity();
            _voronoiMeshSnapshot->importNumber();
            bothDensityAndMass = true;
            break;
    }

    // determine whether to forego the Voronoi mesh
    auto config = find<Configuration>();
    if (bothDensityAndMass && !config->mediaNeedGeneratePosition() && !config->snapshotsNeedGetEntities())
        _voronoiMeshSnapshot->foregoVoronoiMesh();

    // set the domain extent
    _voronoiMeshSnapshot->setExtent(domain());
    return _voronoiMeshSnapshot;
}

////////////////////////////////////////////////////////////////////

VoronoiMeshSnapshot* VoronoiMeshMedium::voronoiMesh() const
{
    return _voronoiMeshSnapshot;
}

////////////////////////////////////////////////////////////////////
