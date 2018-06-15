/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "PointSource.hpp"
#include "PhotonPacket.hpp"
#include "Random.hpp"

//////////////////////////////////////////////////////////////////////

int PointSource::geometryDimension() const
{
    int dimension = 1;
    if (positionZ()) dimension = 2;
    if (positionX() || positionY()) dimension = 3;
    return dimension;
}

//////////////////////////////////////////////////////////////////////

void PointSource::launchNormalized(PhotonPacket* pp, size_t historyIndex, double lambda, double Lw,
                                       RedshiftInterface* rsi) const
{
    // get the source position
    Position bfr(positionX(), positionY(), positionZ());

    // launch the photon packet with isotropic direction
    pp->launch(historyIndex, lambda, Lw, bfr, random()->direction(), rsi);
}

//////////////////////////////////////////////////////////////////////
