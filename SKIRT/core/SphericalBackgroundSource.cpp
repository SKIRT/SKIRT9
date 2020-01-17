/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "SphericalBackgroundSource.hpp"
#include "AngularDistributionInterface.hpp"
#include "PhotonPacket.hpp"
#include "Position.hpp"
#include "Random.hpp"

////////////////////////////////////////////////////////////////////

int SphericalBackgroundSource::intrinsicDimension() const
{
    return 1;
}

////////////////////////////////////////////////////////////////////

namespace
{
    // this function generates a random launch direction for a given launch position on the unit sphere
    // (the launch position on the unit sphere is actually specified as a direction)
    Direction generateDirectionForLaunchPosition(Random* random, Direction bfu)
    {
        // picking a random (theta',phi')
        double thetap = M_PI - acos(sqrt(random->uniform()));
        double phip = 2.0 * M_PI * random->uniform();
        Direction bfkp(thetap, phip);
        double kpx, kpy, kpz;
        bfkp.cartesian(kpx, kpy, kpz);

        // conversion to the regular coordinate system
        double theta, phi;
        bfu.spherical(theta, phi);
        double costheta = cos(theta);
        double sintheta = sin(theta);
        double cosphi = cos(phi);
        double sinphi = sin(phi);
        double kx = costheta * cosphi * kpx - sinphi * kpy + sintheta * cosphi * kpz;
        double ky = costheta * sinphi * kpx + cosphi * kpy + sintheta * sinphi * kpz;
        double kz = -sintheta * kpx + costheta * kpz;
        return Direction(kx, ky, kz);
    }

    // this function returns the normalized probability of launching in a certain direction
    // for a given launch position on the unit sphere (which is actually specified as a direction)
    double probabilityOfDirectionForLaunchPosition(Direction bfk, Direction bfu)
    {
        double costhetap = Vec::dot(bfk, bfu);
        if (costhetap < 0.)
            return -4. * costhetap;
        else
            return 0.;
    }

    // this class serves up the angular distribution interface for a particular launch position on the unit sphere;
    // an instance of this class, configured with a launch position, is handed to each photon packet upon launch
    // so that peel-off photon packets can retrieve the appropriate bias factor for the instrument direction
    class LocalAngularDistribution : public AngularDistributionInterface
    {
    public:
        void setLaunchPosition(Direction bfu) { _bfu = bfu; }
        double probabilityForDirection(Direction bfk) const override
        {
            return probabilityOfDirectionForLaunchPosition(bfk, _bfu);
        }

    private:
        Direction _bfu;
    };

    // setup a local angular distribution instance for each parallel execution thread; this works even if
    // there are multiple sources of this type because each thread handles a single photon packet at a time
    thread_local LocalAngularDistribution _lad;
}

////////////////////////////////////////////////////////////////////

void SphericalBackgroundSource::launchSpecialty(PhotonPacket* pp, size_t historyIndex, double lambda, double Lw,
                                                VelocityInterface* bvi) const
{
    // generate a random intrinsic launch "position" on the unit sphere
    Direction bfu = random()->direction();

    // generate a random launch direction for that position
    Direction bfk = generateDirectionForLaunchPosition(random(), bfu);

    // scale and translate the launch position according to the stellar radius and center
    Position bfr(backgroundRadius(), bfu);
    bfr += Vec(centerX(), centerY(), centerZ());

    // configure the local angular distribution object with the intrinsic launch position
    _lad.setLaunchPosition(bfu);

    // launch the photon packet with the appropriate angular distribution
    pp->launch(historyIndex, lambda, Lw, bfr, bfk, bvi, &_lad);
}

////////////////////////////////////////////////////////////////////
