/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "GeometricSource.hpp"
#include "PhotonPacket.hpp"
#include "Random.hpp"

//////////////////////////////////////////////////////////////////////

int GeometricSource::geometryDimension() const
{
    return geometry()->dimension();
}

//////////////////////////////////////////////////////////////////////

void GeometricSource::launchNormalized(PhotonPacket* pp, size_t historyIndex, double lambda, double Lw,
                                       VelocityInterface* bvi) const
{
    // generate a random position from the geometry
    Position bfr = _geometry->generatePosition();

    // launch the photon packet with isotropic direction
    pp->launch(historyIndex, lambda, Lw, bfr, random()->direction(), bvi);
}

//////////////////////////////////////////////////////////////////////
