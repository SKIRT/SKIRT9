/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "SpecialtySource.hpp"
#include "Configuration.hpp"

//////////////////////////////////////////////////////////////////////

void SpecialtySource::setupSelfBefore()
{
    NormalizedSource::setupSelfBefore();

    // if we have a nonzero bulk velocity, set the interface object to ourselves; otherwise leave it at null pointer
    if (hasVelocity()) _bvi = this;
}

//////////////////////////////////////////////////////////////////////

int SpecialtySource::dimension() const
{
    int velocityDimension = 1;
    if (hasVelocity())
    {
        if (velocityZ()) velocityDimension = 2;
        if (velocityX() || velocityY()) velocityDimension = 3;
    }
    return max(geometryDimension(), velocityDimension);
}

//////////////////////////////////////////////////////////////////////

bool SpecialtySource::hasVelocity() const
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

Vec SpecialtySource::velocity() const
{
    return Vec(velocityX(), velocityY(), velocityZ());
}

//////////////////////////////////////////////////////////////////////

void SpecialtySource::launchNormalized(PhotonPacket* pp, size_t historyIndex, double lambda, double Lw) const
{
    launchSpecialty(pp, historyIndex, lambda, Lw, _bvi);
}

//////////////////////////////////////////////////////////////////////
