/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "GeometricSource.hpp"
#include "Configuration.hpp"
#include "PhotonPacket.hpp"
#include "Random.hpp"
#include "VelocityInterface.hpp"

//////////////////////////////////////////////////////////////////////

void GeometricSource::setupSelfBefore()
{
    NormalizedSource::setupSelfBefore();

    _hasVelocity = hasVelocity();
}

//////////////////////////////////////////////////////////////////////

int GeometricSource::dimension() const
{
    int velocityDimension = hasVelocity() ? velocityDistribution()->dimension() : 1;
    return max(geometry()->dimension(), velocityDimension);
}

//////////////////////////////////////////////////////////////////////

bool GeometricSource::hasVelocity() const
{
    if (velocityDistribution() && velocityMagnitude())
    {
        // refuse velocity for oligochromatic simulations
        // (this function is called from Configure so we cannot precompute this during setup)
        auto config = find<Configuration>();
        if (!config->oligochromatic()) return true;
    }
    return false;
}

//////////////////////////////////////////////////////////////////////

namespace
{
    // an instance of this class offers the velocity interface for a given position
    class LaunchVelocity : public VelocityInterface
    {
    private:
        Vec _bfv;

    public:
        LaunchVelocity() {}
        void setBulkVelocity(Vec bfv) { _bfv = bfv; }
        Vec velocity() const override { return _bfv; }
    };

    // setup a velocity instance (with the redshift interface) for each parallel execution thread; this works even if
    // there are multiple sources of this type because each thread handles a single photon packet at a time
    thread_local LaunchVelocity t_velocity;
}

//////////////////////////////////////////////////////////////////////

void GeometricSource::launchNormalized(PhotonPacket* pp, size_t historyIndex, double lambda, double Lw) const
{
    // generate a random position from the geometry
    Position bfr = _geometry->generatePosition();

    // provide a redshift interface for the appropriate velocity, if applicable
    VelocityInterface* bvi = nullptr;
    if (_hasVelocity)
    {
        Vec bfv = velocityMagnitude() * velocityDistribution()->vector(bfr);
        t_velocity.setBulkVelocity(bfv);
        bvi = &t_velocity;
    }

    // launch the photon packet with isotropic direction
    pp->launch(historyIndex, lambda, Lw, bfr, random()->direction(), bvi);
}

//////////////////////////////////////////////////////////////////////
