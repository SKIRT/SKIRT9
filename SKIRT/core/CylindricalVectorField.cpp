/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "CylindricalVectorField.hpp"

//////////////////////////////////////////////////////////////////////

int CylindricalVectorField::dimension() const
{
    return 3;
}

//////////////////////////////////////////////////////////////////////

Vec CylindricalVectorField::vector(Position bfr) const
{
    // get the radial distance and a unit vector in the rotation direction
    // except that, if at the z-axis, always immediately return a null vector
    Vec u(-bfr.y(), bfr.x(), 0.);
    double r = u.norm();
    if (r == 0) return Vec();
    u /= r;

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
