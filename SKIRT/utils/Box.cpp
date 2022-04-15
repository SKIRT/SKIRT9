/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "Box.hpp"

//////////////////////////////////////////////////////////////////////

bool Box::intersects(Vec r, const Vec k, double& smin, double& smax) const
{
    constexpr double eps = 1e-9;
    smin = -std::numeric_limits<double>::infinity();
    smax = +std::numeric_limits<double>::infinity();

    // --- handle YZ planes ---

    // check for ray parallel to plane
    if (abs(k.x()) < eps)
    {
        // if parallel AND outside box: no intersection possible
        if (r.x() < xmin() || r.x() >= xmax()) return false;
    }
    else
    {
        // find intersection distances and put them in increasing order
        double s1 = (xmin() - r.x()) / k.x();
        double s2 = (xmax() - r.x()) / k.x();
        if (s1 > s2) std::swap(s1, s2);

        // compare with current values
        if (s1 > smin) smin = s1;
        if (s2 < smax) smax = s2;

        // check if ray misses entirely
        if (smin >= smax || smax <= 0) return false;
    }

    // --- handle XZ planes ---

    if (abs(k.y()) < eps)
    {
        if (r.y() < ymin() || r.y() >= ymax()) return false;
    }
    else
    {
        double s1 = (ymin() - r.y()) / k.y();
        double s2 = (ymax() - r.y()) / k.y();
        if (s1 > s2) std::swap(s1, s2);
        if (s1 > smin) smin = s1;
        if (s2 < smax) smax = s2;
        if (smin >= smax || smax <= 0) return false;
    }

    // --- handle XY planes ---

    if (abs(k.z()) < eps)
    {
        if (r.z() < zmin() || r.z() >= zmax()) return false;
    }
    else
    {
        double s1 = (zmin() - r.z()) / k.z();
        double s2 = (zmax() - r.z()) / k.z();
        if (s1 > s2) std::swap(s1, s2);
        if (s1 > smin) smin = s1;
        if (s2 < smax) smax = s2;
        if (smin >= smax || smax <= 0) return false;
    }

    // --- box is definitely intersected ---

    // if origin is inside the box, set first intersection distance to zero
    if (smin < 0.) smin = 0.;
    return true;
}

//////////////////////////////////////////////////////////////////////

namespace
{
    // returns the square of the argument
    double square(double value) { return value * value; }
}

//////////////////////////////////////////////////////////////////////

bool Box::intersects(Vec rc, double r) const
{
    double squaredist = square(r);

    if (rc.x() < xmin())
        squaredist -= square(rc.x() - xmin());
    else if (rc.x() > xmax())
        squaredist -= square(rc.x() - xmax());
    if (rc.y() < ymin())
        squaredist -= square(rc.y() - ymin());
    else if (rc.y() > ymax())
        squaredist -= square(rc.y() - ymax());
    if (rc.z() < zmin())
        squaredist -= square(rc.z() - zmin());
    else if (rc.z() > zmax())
        squaredist -= square(rc.z() - zmax());

    return squaredist >= 0.;
}

//////////////////////////////////////////////////////////////////////
