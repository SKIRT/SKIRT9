/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "CylindricalCell.hpp"
#include "Box.hpp"
#include "Position.hpp"
#include "Quadrics.hpp"
#include <array>

//////////////////////////////////////////////////////////////////////

namespace
{
    // snap a cosine or sine value to 0 or 1 to avoid rounding errors in the calculations
    double snap(double v)
    {
        constexpr double eps = 1e-12;
        if (abs(v) < eps) return 0.;
        if (v > 1. - eps) return 1.;
        if (v < -1. + eps) return -1.;
        return v;
    }
}

//////////////////////////////////////////////////////////////////////

CylindricalCell::CylindricalCell(double Rmin, double phimin, double zmin, double Rmax, double phimax, double zmax)
    : _Rmin(Rmin), _phimin(phimin), _zmin(zmin), _Rmax(Rmax), _phimax(phimax), _zmax(zmax),
      _cosphimin(snap(cos(phimin))), _sinphimin(snap(sin(phimin))), _cosphimax(snap(cos(phimax))),
      _sinphimax(snap(sin(phimax)))
{}

//////////////////////////////////////////////////////////////////////

double CylindricalCell::volume() const
{
    return 0.5 * (_Rmax - _Rmin) * (_Rmax + _Rmin) * (_phimax - _phimin) * (_zmax - _zmin);
}

//////////////////////////////////////////////////////////////////////

Position CylindricalCell::center() const
{
    double R = 0.5 * (_Rmin + _Rmax);
    double phi = 0.5 * (_phimin + _phimax);
    double z = 0.5 * (_zmin + _zmax);
    return Position(R, phi, z, Position::CoordinateSystem::CYLINDRICAL);
}

//////////////////////////////////////////////////////////////////////

bool CylindricalCell::contains(Vec r) const
{
    if (r.z() < _zmin || r.z() >= _zmax) return false;

    double R = sqrt(r.x() * r.x() + r.y() * r.y());
    if (R < _Rmin || R >= _Rmax) return false;

    // as an exception, don't exclude pi from the maximum
    double phi = atan2(r.y(), r.x());
    if (phi < _phimin || phi > _phimax) return false;
    if (phi != M_PI && phi == _phimax) return false;

    return true;
}

//////////////////////////////////////////////////////////////////////

Box CylindricalCell::boundingBox() const
{
    // the (x,y) coordinates of the four corner points
    double x1 = _Rmin * _cosphimin;
    double x2 = _Rmin * _cosphimax;
    double x3 = _Rmax * _cosphimin;
    double x4 = _Rmax * _cosphimax;
    double y1 = _Rmin * _sinphimin;
    double y2 = _Rmin * _sinphimax;
    double y3 = _Rmax * _sinphimin;
    double y4 = _Rmax * _sinphimax;

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

double CylindricalCell::intersection(Vec r, const Vec k) const
{
    // allocate room for the 8 possible intersections (1 per plane and 2 per cylinder)
    // plus the starting position (which starts the first segment).
    // initialize the array to zeroes:
    //  - leave the first value at zero to represent the starting position
    //  - overwrite the other values for each intersection, or leave at zero if there is none
    enum { START, ZMIN, ZMAX, PHIMIN, PHIMAX, RMIN1, RMIN2, RMAX1, RMAX2, LEN };
    std::array<double, LEN> sv = {};

    // intersections with horizontal planes
    sv[ZMIN] = Quadrics::intersectionHorizontalPlane(r, k, zmin());
    sv[ZMAX] = Quadrics::intersectionHorizontalPlane(r, k, zmax());

    // intersections with meridional planes
    sv[PHIMIN] = Quadrics::intersectionMeridionalPlane(r, k, _sinphimin, _cosphimin);
    sv[PHIMAX] = Quadrics::intersectionMeridionalPlane(r, k, _sinphimax, _cosphimax);

    // intersections with boundary cylinders
    double kq2 = k.x() * k.x() + k.y() * k.y();
    Quadrics::distinctIntersectionsCylinder(r, k, kq2, _Rmin, sv[RMIN1], sv[RMIN2]);
    Quadrics::distinctIntersectionsCylinder(r, k, kq2, _Rmax, sv[RMAX1], sv[RMAX2]);

    // sort the intersection points
    std::sort(sv.begin(), sv.end());

    // accumulate the length of all segments that are beyond the starting position and inside the cell
    // (there is always at least one zero in the list)
    double length = 0.;
    for (auto sp = std::find_if(sv.begin() + 1, sv.end(), [](double s) { return s > 0; }); sp != sv.end(); ++sp)
    {
        // get the start and end points for this segment
        double s1 = *(sp - 1);
        double s2 = *sp;

        // if the midpoint is inside the cell, add the segment length
        double s = 0.5 * (s1 + s2);
        if (contains(Vec(r.x() + s * k.x(), r.y() + s * k.y(), r.z() + s * k.z()))) length += s2 - s1;
    }
    return length;
}

//////////////////////////////////////////////////////////////////////
