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
    _adaptiveMeshSnapshot = new AdaptiveMeshSnapshot;
    _adaptiveMeshSnapshot->open(this, filename(), "adaptive mesh cells");

    // honor custom column reordering
    _adaptiveMeshSnapshot->useColumns(useColumns());

    // configure the mass or density column
    switch (massType())
    {
        case MassType::MassDensity: _adaptiveMeshSnapshot->importMassDensity(); break;
        case MassType::Mass: _adaptiveMeshSnapshot->importMass(); break;
        case MassType::MassDensityAndMass:
            _adaptiveMeshSnapshot->importMassDensity();
            _adaptiveMeshSnapshot->importMass();
            break;

        case MassType::NumberDensity: _adaptiveMeshSnapshot->importNumberDensity(); break;
        case MassType::Number: _adaptiveMeshSnapshot->importNumber(); break;
        case MassType::NumberDensityAndNumber:
            _adaptiveMeshSnapshot->importNumberDensity();
            _adaptiveMeshSnapshot->importNumber();
            break;
    }

    // set the domain extent
    _adaptiveMeshSnapshot->setExtent(domain());
    return _adaptiveMeshSnapshot;
}

////////////////////////////////////////////////////////////////////

AdaptiveMeshSnapshot* AdaptiveMeshMedium::adaptiveMesh() const
{
    return _adaptiveMeshSnapshot;
}

////////////////////////////////////////////////////////////////////
