/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef CONVEXPOLYGON_HPP
#define CONVEXPOLYGON_HPP

#include "Basics.hpp"
#include <array>

////////////////////////////////////////////////////////////////////

/** The ConvexPolygon class represents a planar convex polygon, and provides support for extracting
    specific properties such as the bottom-most and top-most points. There is a (fairly small)
    compile-time limit on the number of points to avoid memory allocation on the heap. The user
    should ensure that at least 3 and at most \em maxPoints points are added to the polygon,
    that they are added in consecutive clockwise order, that there is a nontrivial distance between
    all points, and that the resulting polygon is convex. If these restrictions are violated, the
    behavior of some of the functions in this class is undefined. */
class ConvexPolygon
{
    // ================== Constants ==================
public:
    /** The maximum number of points that can be held by a polygon object. */
    constexpr static size_t maxPoints = 8;

    // ================== Constructing ==================

public:
    /** Constructs an empty polygon. Use the add() function to add points to the polygon. */
    ConvexPolygon();

    /** Adds a point to the polygon. The user should ensure that at least 3 and at most 8 points
        are added to the polygon, that they are added in consecutive clockwise order, that there is
        a nontrivial distance between all points, and that the resulting polygon is convex. */
    void add(double x, double y);

    // ================== Getting points and properties ==================

public:
    /** Returns the number of points in the polygon. */
    size_t count() const { return _n; }

    /** Returns the x-coordinate of the point at the specified zero-based index. */
    double x(size_t i) const { return _xv[i]; }

    /** Returns the y-coordinate of the point at the specified zero-based index. */
    double y(size_t i) const { return _yv[i]; }

    /** Returns the top-most y-coordinate of the polygon. */
    double top() const { return _yv[_topIndex]; }

    /** Returns the bottom-most y-coordinate of the polygon. */
    double bottom() const { return _yv[_botIndex]; }

    /** Returns the x-coordinate for the left-hand intersection point of the polygon with the
        horizontal line at the specified y-coordinate. If the specified y-coordinate is out
        of range, the x-coordinate of the top or bottom point is returned. */
    double leftFor(double y) const;

    /** Returns the x-coordinate for the right-hand intersection point of the polygon with the
        horizontal line at the specified y-coordinate. If the specified y-coordinate is out
        of range, the x-coordinate of the top or bottom point is returned. */
    double rightFor(double y) const;

    // ================== Data members ==================

private:
    std::array<double, maxPoints> _xv;  // the x coordinates of the points in the polygon
    std::array<double, maxPoints> _yv;  // the y coordinates of the points in the polygon
    size_t _n = 0;                      // the number of points in the polygon
    size_t _topIndex = 0;               // the index of the top-most point (right-most if there is more than one)
    size_t _botIndex = 0;               // the index of the bottom-most point (left-most if there is more than one)
};

////////////////////////////////////////////////////////////////////

#endif
