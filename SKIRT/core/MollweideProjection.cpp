/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "MollweideProjection.hpp"

////////////////////////////////////////////////////////////////////

void MollweideProjection::fromGlobeToRectangle(double longitude, double latitude, double& x, double& y) const
{
    const double eps = M_PI / 500;  // arbitrary accuracy choice

    double theta = latitude;
    if (theta > -M_PI_2 + eps && theta < M_PI_2 - eps)
    {
        while (true)
        {
            double diff = (2 * theta + sin(2 * theta) - M_PI * sin(latitude)) / (2 + 2 * cos(2 * theta));
            theta -= diff;
            if (abs(diff) < eps) break;
        }
    }
    x = longitude / M_PI * cos(theta);
    y = sin(theta);
}

////////////////////////////////////////////////////////////////////

bool MollweideProjection::fromRectangleToGlobe(double x, double y, double& longitude, double& latitude) const
{
    double alpha = asin(y);
    longitude = M_PI * x / cos(alpha);
    if (longitude >= -M_PI && longitude <= M_PI)
    {
        latitude = asin((2 * alpha + sin(2 * alpha)) / M_PI);
        return true;
    }
    return false;
}

////////////////////////////////////////////////////////////////////
