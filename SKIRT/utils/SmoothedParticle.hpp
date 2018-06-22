/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef SMOOTHEDPARTICLE_HPP
#define SMOOTHEDPARTICLE_HPP

#include "Vec.hpp"
class Box;

////////////////////////////////////////////////////////////////////

/** SmoothedParticle is an almost trivial helper class for working with smoothed particles, usually
    imported from a smoothed particle hydrodynamical (SPH) simulation. A SmoothedParticle object
    holds the particle's center position, smoothing length, and effective mass (after any
    multipliers have been applied). Isolating these properties in a simple object allows efficient
    storage and retrieval, for example when organizing smoothed particles in a search grid. */
class SmoothedParticle
{
public:
    /** The constructor arguments specify the particle attributes: coordinates of the center,
        smoothing length, and effective mass. */
    SmoothedParticle(double x, double y, double z, double h, double M)
        : _r{x,y,z}, _h(h), _M(M) { }

    /** This function returns the coordinates of the center of the particle. */
    Vec center() const { return Vec(_r[0], _r[1], _r[2]); }

    /** This function returns the x, y, or z-coordinate of the center of the particle, depending on
        the value of \em dir: 1->x, 2->y, 3->z. For any other value of \em dir the behavior is
        undefined. */
    double center(int dir) const { return _r[dir-1]; }

    /** This function returns the smoothing length of the particle. */
    double radius() const { return _h; }

    /** This function returns the mass of the particle. */
    double mass() const { return _M; }

private:
    // data members received as constructor arguments
    double _r[3];   // center coordinates
    double _h;      // smoothing length
    double _M;      // total mass
};

////////////////////////////////////////////////////////////////////

#endif
