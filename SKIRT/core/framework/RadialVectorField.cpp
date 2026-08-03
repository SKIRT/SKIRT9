/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "RadialVectorField.hpp"

//////////////////////////////////////////////////////////////////////

int RadialVectorField::dimension() const
{
    return 3;
}

//////////////////////////////////////////////////////////////////////

Vec RadialVectorField::vector(Position bfr) const
{
    // get the radial distance and a unit vector in the radial direction
    // except that, if at the origin, always immediately return a null vector
    double r = bfr.norm();
    if (r == 0) return Vec();
    Vec u = bfr / r;

    // set the magnitude to a default value of unity and update if needed
    double v = 1.;
    if (unityRadius() > 0.)
    {
        if ((exponent() > 0. && r < unityRadius()) || (exponent() < 0. && r > unityRadius()))
            v = pow(r / unityRadius(), exponent());
    }

    // return a vector with proper magnitude and direction
    return v * u;
}

//////////////////////////////////////////////////////////////////////
