/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "ConvexPolygon.hpp"

////////////////////////////////////////////////////////////////////

ConvexPolygon::ConvexPolygon() {}

////////////////////////////////////////////////////////////////////

void ConvexPolygon::add(double x, double y)
{
    if (_n < maxPoints)
    {
        _xv[_n] = x;
        _yv[_n] = y;
        if (y >= _yv[_topIndex]) _topIndex = _n;
        if (y <= _yv[_botIndex]) _botIndex = _n;
        ++_n;
    }
}

////////////////////////////////////////////////////////////////////

namespace
{
    // returns the x-coordinate of the intersection of the horizontal line at the specified y-coordinate
    // with the line through the specified points p1 and p2
    double intersect(double y, double x1, double y1, double x2, double y2)
    {
        return x1 + (x2 - x1) * (y - y1) / (y2 - y1);
    }
}

double ConvexPolygon::leftFor(double y) const
{
    // handle out-of-range cases
    if (y > _yv[_topIndex]) return _xv[_topIndex];
    if (y <= _yv[_botIndex]) return _xv[_botIndex];

    // loop through points in clockwise order starting at bottom point
    // and identify the point ending the intersection segment
    size_t i = _botIndex + 1;
    while (y > _yv[i % _n]) ++i;

    // return the intersection point with that segment
    size_t i1 = (i - 1) % _n;
    size_t i2 = i % _n;
    return intersect(y, _xv[i1], _yv[i1], _xv[i2], _yv[i2]);
}

////////////////////////////////////////////////////////////////////

double ConvexPolygon::rightFor(double y) const
{
    // handle out-of-range cases
    if (y < _yv[_botIndex]) return _xv[_botIndex];
    if (y >= _yv[_topIndex]) return _xv[_topIndex];

    // loop through points in clockwise order starting at top point
    // and identify the point ending the intersection segment
    size_t i = _topIndex + 1;
    while (y < _yv[i % _n]) ++i;

    // return the intersection point with that segment
    size_t i1 = (i - 1) % _n;
    size_t i2 = i % _n;
    return intersect(y, _xv[i1], _yv[i1], _xv[i2], _yv[i2]);
}

////////////////////////////////////////////////////////////////////
