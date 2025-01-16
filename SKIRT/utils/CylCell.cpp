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

    // the cosine and sine for each of the azimuthal angles
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

namespace
{
    // determines the solutions of x^2 + 2*b*x + c = 0
    // - if there are two distinct real solutions, they are stored in x1,x2 and the function returns true
    // - otherwise, x1 and x2 remain unchanged and the function returns false
    bool solutions(double b, double c, double& x1, double& x2)
    {
        // x1 == -b - sqrt(b*b-c)
        // x2 == -b + sqrt(b*b-c)
        // x1*x2 == c

        if (b * b > c)  // if discriminant is strictly positive, there are two distinct real solutions
        {
            if (b > 0)  // x1 is always negative
            {
                x1 = -b - sqrt(b * b - c);
                x2 = c / x1;
            }
            else  // x2 is always positive
            {
                x2 = -b + sqrt(b * b - c);
                x1 = c / x2;
            }
            return true;
        }
        return false;
    }
}

//////////////////////////////////////////////////////////////////////

double CylCell::intersection(Vec r, const Vec k) const
{
    constexpr double eps = 1e-9;
    double smin = -std::numeric_limits<double>::infinity();
    double smax = +std::numeric_limits<double>::infinity();

    // --- handle horizontal planes ---

    // check for ray parallel to plane
    if (abs(k.z()) < eps)
    {
        // if parallel AND outside: no intersection possible
        if (r.z() < zmin() || r.z() >= zmax()) return 0.;
    }
    else
    {
        // find intersection distances and put them in increasing order
        smin = (zmin() - r.z()) / k.z();
        smax = (zmax() - r.z()) / k.z();
        if (smin > smax) std::swap(smin, smax);

        // check if ray misses entirely
        if (smax <= 0) return 0.;
    }

    // --- handle meridional phimin plane ---
    {
        double cosphi = cos(_phimin);
        double sinphi = sin(_phimin);
        double q = k.x() * sinphi - k.y() * cosphi;
        if (abs(q) >= eps)
        {
            double s = (r.x() * sinphi - r.y() * cosphi) / q;
            smin = min(smin, s);
            smax = max(smax, s);
        }
    }

    // --- handle meridional phimax plane ---
    {
        double cosphi = cos(_phimax);
        double sinphi = sin(_phimax);
        double q = k.x() * sinphi - k.y() * cosphi;
        if (abs(q) >= eps)
        {
            double s = (r.x() * sinphi - r.y() * cosphi) / q;
            smin = min(smin, s);
            smax = max(smax, s);
        }
    }

    // --- handle cylinders ---
    {
        double a = k.x() * k.x() + k.y() * k.y();
        if (abs(a) >= eps)
        {
            // outer
            double b = (r.x() * k.x() + r.y() * k.y()) / a;
            double c = (r.x() * r.x() + r.y() * r.y() - _Rmax * _Rmax) / a;
            double s1, s2;
            if (solutions(b, c, s1, s2))
            {
                smin = min({smin, s1, s2});
                smax = max({smax, s1, s2});
            }

            // inner
            c = (r.x() * r.x() + r.y() * r.y() - _Rmin * _Rmin) / a;
            if (solutions(b, c, s1, s2))
            {
                smin = min({smin, s1, s2});
                smax = max({smax, s1, s2});
            }
        }
    }

    // TODO: handles cases where the line lies inside a border plane or cylinder
    // TODO: correct for segment outside of inner cylinder, if applicable

    // if origin is inside the cell, set first intersection distance to zero
    if (smin < 0.) smin = 0.;
    return smax - smin;
}

//////////////////////////////////////////////////////////////////////
