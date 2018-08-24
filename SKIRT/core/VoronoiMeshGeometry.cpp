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
    auto snapshot = new VoronoiMeshSnapshot;
    snapshot->open(this, filename(), "Voronoi sites");

    // configure the mass or density column (position columns are configured by the snapshot itself)
    switch (massType())
    {
    case MassType::MassDensity: snapshot->importMassDensity(); break;
    case MassType::Mass: snapshot->importMass(); break;
    case MassType::NumberDensity: snapshot->importNumberDensity(); break;
    case MassType::Number: snapshot->importNumber(); break;
    }

    // set the domain extent
    snapshot->setExtent(domain());
    return snapshot;
}

////////////////////////////////////////////////////////////////////
