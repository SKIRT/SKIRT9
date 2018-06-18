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
    // point source symmetry depends on

    // ... the point's location
    int positionDimension = 1;
    if (positionZ()) positionDimension = 2;
    if (positionX() || positionY()) positionDimension = 3;

    // ... the angular distribution of the emission
    int angularDimension = angularDistribution() ? angularDistribution()->dimension() : 1;

    return max(positionDimension, angularDimension);
}

//////////////////////////////////////////////////////////////////////

void PointSource::launchNormalized(PhotonPacket* pp, size_t historyIndex, double lambda, double Lw,
                                       RedshiftInterface* rsi) const
{
    // get the source position
    Position bfr(positionX(), positionY(), positionZ());

    // launch the photon packet with the appropriate angular distribution
    pp->launch(historyIndex, lambda, Lw, bfr,
               angularDistribution() ? angularDistribution()->generateDirection() : random()->direction(),
               rsi, angularDistribution());
}

//////////////////////////////////////////////////////////////////////
