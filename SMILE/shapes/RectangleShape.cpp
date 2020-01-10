/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "RectangleShape.hpp"

////////////////////////////////////////////////////////////////////

void RectangleShape::paintSelf()
{
    double x1 = x() - 0.5 * width();
    double x2 = x() + 0.5 * width();
    double y1 = y() - 0.5 * height();
    double y2 = y() + 0.5 * height();

    drawLine(x1, y1, x1, y2);
    drawLine(x1, y2, x2, y2);
    drawLine(x2, y2, x2, y1);
    drawLine(x2, y1, x1, y1);
}

////////////////////////////////////////////////////////////////////
