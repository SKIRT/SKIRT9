/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "AdaptiveMeshMedium.hpp"
#include "AdaptiveMeshSnapshot.hpp"

////////////////////////////////////////////////////////////////////

Snapshot* AdaptiveMeshMedium::createAndOpenSnapshot()
{
    // create and open the snapshot
    auto snapshot = new AdaptiveMeshSnapshot;
    snapshot->open(this, filename(), "adaptive mesh cells");

    // configure the mass or density column
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
