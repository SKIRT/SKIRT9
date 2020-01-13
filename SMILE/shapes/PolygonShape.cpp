/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "PolygonShape.hpp"

////////////////////////////////////////////////////////////////////

void PolygonShape::paintSelf()
{
    double x1 = radius();
    double y1 = 0.;

    for (int i = 1; i <= numSides(); ++i)
    {
        double theta = 2 * M_PI * i / numSides();
        double x2 = radius() * cos(theta);
        double y2 = radius() * sin(theta);

        if (align())
            drawLine(x() + y1, y() + x1, x() + y2, y() + x2);
        else
            drawLine(x() + x1, y() + y1, x() + x2, y() + y2);

        x1 = x2;
        y1 = y2;
    }
}

////////////////////////////////////////////////////////////////////
