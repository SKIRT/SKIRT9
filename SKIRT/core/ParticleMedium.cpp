/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "ParticleMedium.hpp"
#include "ParticleSnapshot.hpp"

////////////////////////////////////////////////////////////////////

Snapshot* ParticleMedium::createAndOpenSnapshot()
{
    // create and open the snapshot
     _particleSnapshot = new ParticleSnapshot;
    _particleSnapshot->open(this, filename(), "smoothed particles");

    // honor custom column reordering
    _particleSnapshot->useColumns(useColumns());

    // configure the position and size columns
    _particleSnapshot->importPosition();
    _particleSnapshot->importSize();

    // configure the mass or number column
    switch (massType())
    {
        case MassType::Mass: _particleSnapshot->importMass(); break;
        case MassType::Number: _particleSnapshot->importNumber(); break;
    }

    // set the smoothing kernel
    _particleSnapshot->setSmoothingKernel(smoothingKernel());

    return _particleSnapshot;
}

////////////////////////////////////////////////////////////////////

double ParticleMedium::massInBox(const Box& box) const
{
    return _particleSnapshot->massInBox(box);
}


////////////////////////////////////////////////////////////////////
