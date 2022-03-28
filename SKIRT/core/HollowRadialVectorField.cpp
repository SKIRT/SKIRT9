/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "HollowRadialVectorField.hpp"

//////////////////////////////////////////////////////////////////////

int HollowRadialVectorField::dimension() const
{
    return 3;
}

//////////////////////////////////////////////////////////////////////

Vec HollowRadialVectorField::vector(Position bfr) const
{
    // get the radial distance; if below r_0, return a null vector
    double r = bfr.norm();
    if (r <= zeroRadius()) return Vec();

    // get a unit vector in the direction of the current position
    Vec u = bfr / r;

    // calculate the magnitude as a function of radial distance
    double v = pow(1. - zeroRadius() / r, exponent());

    // return a vector with proper magnitude and direction
    return v * u;
}

//////////////////////////////////////////////////////////////////////
