/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "TetraMeshMedium.hpp"
#include "Configuration.hpp"
#include "TetraMeshSnapshot.hpp"

////////////////////////////////////////////////////////////////////

Snapshot* TetraMeshMedium::createAndOpenSnapshot()
{
    // create and open the snapshot
    _tetraMeshSnapshot = new TetraMeshSnapshot;
    _tetraMeshSnapshot->open(this, filename(), "Tetra sites");

    // honor custom column reordering
    _tetraMeshSnapshot->useColumns(useColumns());

    // configure the position columns
    _tetraMeshSnapshot->importPosition();

    // configure the density and/or mass column(s)
    bool bothDensityAndMass = false;
    switch (massType())
    {
        case MassType::MassDensity: _tetraMeshSnapshot->importMassDensity(); break;
        case MassType::Mass: _tetraMeshSnapshot->importMass(); break;
        case MassType::MassDensityAndMass:
            _tetraMeshSnapshot->importMassDensity();
            _tetraMeshSnapshot->importMass();
            bothDensityAndMass = true;
            break;

        case MassType::NumberDensity: _tetraMeshSnapshot->importNumberDensity(); break;
        case MassType::Number: _tetraMeshSnapshot->importNumber(); break;
        case MassType::NumberDensityAndNumber:
            _tetraMeshSnapshot->importNumberDensity();
            _tetraMeshSnapshot->importNumber();
            bothDensityAndMass = true;
            break;
    }

    // determine whether to forego the Tetra mesh
    // auto config = find<Configuration>();
    // if (bothDensityAndMass && !config->mediaNeedGeneratePosition() && !config->snapshotsNeedGetEntities())
    //     _tetraMeshSnapshot->foregoTetraMesh();

    // set the domain extent
    _tetraMeshSnapshot->setExtent(domain());
    return _tetraMeshSnapshot;
}

////////////////////////////////////////////////////////////////////

TetraMeshSnapshot* TetraMeshMedium::tetraMesh() const
{
    return _tetraMeshSnapshot;
}

////////////////////////////////////////////////////////////////////
