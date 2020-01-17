/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "CubicalBackgroundSource.hpp"
#include "AngularDistributionInterface.hpp"
#include "PhotonPacket.hpp"
#include "Position.hpp"
#include "Random.hpp"

////////////////////////////////////////////////////////////////////

int CubicalBackgroundSource::intrinsicDimension() const
{
    return 3;
}

////////////////////////////////////////////////////////////////////

namespace
{
    // this function generates a random launch wall of the cube, specified by its outward normal
    Direction generateLaunchWall(Random* random)
    {
        double X = random->uniform();
        if (X < 1. / 6.) return Direction(-1., 0., 0.);
        if (X < 2. / 6.) return Direction(1., 0., 0.);
        if (X < 3. / 6.) return Direction(0., -1., 0.);
        if (X < 4. / 6.) return Direction(0., 1., 0.);
        if (X < 5. / 6.) return Direction(0., 0., -1.);
        return Direction(0., 0., 1.);
    }

    // this function generates a random launch position on the wall specified by its outward normal
    Position generateLaunchPositionOnLaunchWall(Random* random, Direction bfu)
    {
        double t1 = 2. * random->uniform() - 1.;
        double t2 = 2. * random->uniform() - 1.;
        if (bfu.x()) return Position(bfu.x(), t1, t2);
        if (bfu.y()) return Position(t1, bfu.y(), t2);
        return Position(t1, t2, bfu.z());
    }

    // this function generates a random launch direction for a given launch wall
    // (which is specified by its outward normal)
    Direction generateDirectionForLaunchWall(Random* random, Direction bfu)
    {
        // picking a random (theta',phi')
        double thetap = M_PI - acos(sqrt(random->uniform()));
        double phip = 2.0 * M_PI * random->uniform();
        Direction bfkp(thetap, phip);
        double kpx, kpy, kpz;
        bfkp.cartesian(kpx, kpy, kpz);

        // conversion to the regular coordinate system
        if (bfu.x() == -1.) return Direction(-kpz, -kpy, -kpx);
        if (bfu.x() == 1.) return Direction(kpz, kpy, -kpx);
        if (bfu.y() == -1.) return Direction(kpy, -kpz, -kpx);
        if (bfu.y() == 1.) return Direction(-kpy, kpz, -kpx);
        if (bfu.z() == -1.) return Direction(-kpx, kpy, -kpz);
        return Direction(kpx, kpy, kpz);
    }

    // this function returns the normalized probability of launching in a certain direction
    // for a given launch wall (which is specified by its outward normal)
    double probabilityOfDirectionForLaunchWall(Direction bfk, Direction bfu)
    {
        double costhetap = Vec::dot(bfk, bfu);
        if (costhetap < 0.)
            return -4. * costhetap;
        else
            return 0.;
    }

    // this class serves up the angular distribution interface for a particular launch wall;
    // an instance of this class, configured with a launch wall, is handed to each photon packet upon launch
    // so that peel-off photon packets can retrieve the appropriate bias factor for the instrument direction
    class LocalAngularDistribution : public AngularDistributionInterface
    {
    public:
        void setLaunchWall(Direction bfu) { _bfu = bfu; }
        double probabilityForDirection(Direction bfk) const override
        {
            return probabilityOfDirectionForLaunchWall(bfk, _bfu);
        }

    private:
        Direction _bfu;
    };

    // setup a local angular distribution instance for each parallel execution thread; this works even if
    // there are multiple sources of this type because each thread handles a single photon packet at a time
    thread_local LocalAngularDistribution _lad;
}

////////////////////////////////////////////////////////////////////

void CubicalBackgroundSource::launchSpecialty(PhotonPacket* pp, size_t historyIndex, double lambda, double Lw,
                                              VelocityInterface* bvi) const
{
    // generate a random launch wall
    Direction bfu = generateLaunchWall(random());

    // generate a random launch direction for that position
    Direction bfk = generateDirectionForLaunchWall(random(), bfu);

    // generate a launch position on the selected wall,
    // and scale and translate the position according to the cube's size and center
    Position bfr = generateLaunchPositionOnLaunchWall(random(), bfu);
    bfr *= 0.5 * edgeLength();
    bfr += Vec(centerX(), centerY(), centerZ());

    // configure the local angular distribution object with the intrinsic launch position
    _lad.setLaunchWall(bfu);

    // launch the photon packet with the appropriate angular distribution
    pp->launch(historyIndex, lambda, Lw, bfr, bfk, bvi, &_lad);
}

////////////////////////////////////////////////////////////////////
