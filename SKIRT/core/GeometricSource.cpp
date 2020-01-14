/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "GeometricSource.hpp"
#include "PhotonPacket.hpp"
#include "Configuration.hpp"
#include "Random.hpp"

//////////////////////////////////////////////////////////////////////

void GeometricSource::setupSelfBefore()
{
    NormalizedSource::setupSelfBefore();

    // if we have a nonzero bulk velocity, set the interface object to ourselves; otherwise leave it at null pointer
    if (hasVelocity()) _bvi = this;
}

//////////////////////////////////////////////////////////////////////

int GeometricSource::dimension() const
{
    int velocityDimension = 1;
    if (hasVelocity())
    {
        if (velocityZ()) velocityDimension = 2;
        if (velocityX() || velocityY()) velocityDimension = 3;
    }
    return max(geometry()->dimension(), velocityDimension);
}

//////////////////////////////////////////////////////////////////////

bool GeometricSource::hasVelocity() const
{
    if (velocityX() || velocityY() || velocityZ())
    {
        // refuse velocity for oligochromatic simulations
        // (this function is called from Configure so we cannot precompute this during setup)
        auto config = find<Configuration>();
        if (!config->oligochromatic()) return true;
    }
    return false;
}

//////////////////////////////////////////////////////////////////////

Vec GeometricSource::velocity() const
{
    return Vec(velocityX(), velocityY(), velocityZ());
}

//////////////////////////////////////////////////////////////////////

void GeometricSource::launchNormalized(PhotonPacket* pp, size_t historyIndex, double lambda, double Lw) const
{
    // generate a random position from the geometry
    Position bfr = _geometry->generatePosition();

    // launch the photon packet with isotropic direction
    pp->launch(historyIndex, lambda, Lw, bfr, random()->direction(), _bvi);
}

//////////////////////////////////////////////////////////////////////
