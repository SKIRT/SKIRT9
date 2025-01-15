/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "CylCell.hpp"

//////////////////////////////////////////////////////////////////////

bool CylCell::contains(Vec r) const
{
    double R = sqrt(r.x() * r.x() + r.y() * r.y());
    double phi = atan2(r.y(), r.x());
    double z = r.z();
    return contains(R, phi, z);
}

//////////////////////////////////////////////////////////////////////

Box CylCell::boundingBox() const
{
    // The bounds along the z-axis are the same for Cylindrical and Cartesian coordinates,
    // so we just need to determine the bounding rectangle projected on the xy plane.
    // This bounding rectangle must of course enclose the four corner points of the cell.
    // In addition, if the cell straddles one of the coordinate axes, the bounding rectangle
    // must also enclose a point on that axis at radius Rmax.

    // the cosine and sine for each of the aximuthal bounding planes
    double cosphimin = cos(_phimin);
    double sinphimin = sin(_phimin);
    double cosphimax = cos(_phimax);
    double sinphimax = cos(_phimax);

    // the (x,y) coordinates of the four corner points
    double x1 = _Rmin * cosphimin;
    double x2 = _Rmin * cosphimax;
    double x3 = _Rmax * cosphimin;
    double x4 = _Rmax * cosphimax;
    double y1 = _Rmin * sinphimin;
    double y2 = _Rmin * sinphimax;
    double y3 = _Rmax * sinphimin;
    double y4 = _Rmax * sinphimax;

    // angles for coordinate axis directions (cannot straddle negative x-axis)
    constexpr double negy = -M_PI_2;
    constexpr double posx = 0.;
    constexpr double posy = M_PI_2;

    // min/max coordinates
    double xmin = min({x1, x2, x3, x4});
    double ymin = (_phimin <= negy && _phimax >= negy) ? -_Rmax : min({y1, y2, y3, y4});
    double xmax = (_phimin <= posx && _phimax >= posx) ? _Rmax : max({x1, x2, x3, x4});
    double ymax = (_phimin <= posy && _phimax >= posy) ? _Rmax : max({y1, y2, y3, y4});

    // bounding box
    return Box(xmin, ymin, _zmin, xmax, ymax, _zmax);
}

//////////////////////////////////////////////////////////////////////

double CylCell::intersection(Vec r, const Vec k) const
{
    constexpr double eps = 1e-9;
    double smin = -std::numeric_limits<double>::infinity();
    double smax = +std::numeric_limits<double>::infinity();

    // --- handle XY planes ---

    // check for ray parallel to plane
    if (abs(k.z()) < eps)
    {
        // if parallel AND outside: no intersection possible
        if (r.z() < zmin() || r.z() >= zmax()) return 0.;
    }
    else
    {
        // find intersection distances and put them in increasing order
        double s1 = (zmin() - r.z()) / k.z();
        double s2 = (zmax() - r.z()) / k.z();
        if (s1 > s2) std::swap(s1, s2);

        // compare with current values
        if (s1 > smin) smin = s1;
        if (s2 < smax) smax = s2;

        // check if ray misses entirely
        if (smin >= smax || smax <= 0) return 0.;
    }

    return 0.;
}

//////////////////////////////////////////////////////////////////////
