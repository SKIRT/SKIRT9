/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "HammerAitoffProjection.hpp"

////////////////////////////////////////////////////////////////////

void HammerAitoffProjection::fromGlobeToRectangle(double longitude, double latitude, double& x, double& y) const
{
    double t = 1 / sqrt(1 + cos(latitude) * cos(longitude / 2));
    x = t * cos(latitude) * sin(longitude / 2);
    y = t * sin(latitude);
}

////////////////////////////////////////////////////////////////////

bool HammerAitoffProjection::fromRectangleToGlobe(double x, double y, double& longitude, double& latitude) const
{
    double r2 = x * x + y * y;
    if (r2 < 1)
    {
        double q = sqrt(2 - r2);
        longitude = 2 * atan(x * q / (q * q - 1));
        latitude = asin(y * q);
        return true;
    }
    return false;
}

////////////////////////////////////////////////////////////////////
