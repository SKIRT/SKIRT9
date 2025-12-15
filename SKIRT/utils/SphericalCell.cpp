/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "SphericalCell.hpp"
#include "Box.hpp"
#include "Cubic.hpp"
#include "Position.hpp"
#include "Quadratic.hpp"
#include <array>

//////////////////////////////////////////////////////////////////////

SphericalCell::SphericalCell(double rmin, double thetamin, double phimin, double rmax, double thetamax, double phimax)
    : _rmin(rmin), _thetamin(thetamin), _phimin(phimin), _rmax(rmax), _thetamax(thetamax), _phimax(phimax),
      _cosphimin(cos(phimin)), _sinphimin(sin(phimin)), _cosphimax(cos(phimax)), _sinphimax(sin(phimax)),
      _costhetamin(cos(thetamin)), _sinthetamin(sin(thetamin)), _costhetamax(cos(thetamax)), _sinthetamax(sin(thetamax))
{}

//////////////////////////////////////////////////////////////////////

double SphericalCell::volume() const
{
    return (1. / 3.) * Cubic::pow3(_rmin, _rmax) * (_costhetamin - _costhetamax) * (_phimax - _phimin);
}

//////////////////////////////////////////////////////////////////////

Position SphericalCell::center() const
{
    double r = 0.5 * (_rmin + _rmax);
    double theta = 0.5 * (_thetamin + _thetamax);
    double phi = 0.5 * (_phimin + _phimax);
    return Position(r, theta, phi, Position::CoordinateSystem::SPHERICAL);
}

//////////////////////////////////////////////////////////////////////

bool SphericalCell::contains(Vec bfr) const
{
    double r = bfr.norm();
    if (r < _rmin || r >= _rmax) return false;
    // avoid trying to determine the angles for the origin
    if (r == 0. && _rmin == 0.) return true;

    // as an exception, don't exclude pi from the maximum
    double theta = acos(bfr.z() / r);
    if (theta < _thetamin || theta > _thetamax) return false;
    if (theta != M_PI && theta == _thetamax) return false;

    // as an exception, don't exclude pi from the maximum
    double phi = atan2(bfr.y(), bfr.x());
    if (phi < _phimin || phi > _phimax) return false;
    if (phi != M_PI && phi == _phimax) return false;

    return true;
}

//////////////////////////////////////////////////////////////////////

Box SphericalCell::boundingBox() const
{
    // find the vertical bounds in meridional plane
    double z1 = _rmin * _costhetamin;
    double z2 = _rmin * _costhetamax;
    double z3 = _rmax * _costhetamin;
    double z4 = _rmax * _costhetamax;
    double zmin = min({z1, z2, z3, z4});
    double zmax = max({z1, z2, z3, z4});

    // trivial class representing an equatorial Cartesian bounding box enclosing the projection of points defined
    // in spherical coordinates; because the sine and cosine for the relevant angles are precalculated, these
    // values are specified rather than the angles themselves.
    struct EqBox
    {
        double xmin, ymin, xmax, ymax;

        EqBox(double r, double sintheta, double cosphi, double sinphi)
            : xmin(r * sintheta * cosphi), ymin(r * sintheta * sinphi), xmax(xmin), ymax(ymin)
        {}

        void extend(double r, double sintheta, double cosphi, double sinphi)
        {
            double x = r * sintheta * cosphi;
            double y = r * sintheta * sinphi;
            xmin = min(xmin, x);
            ymin = min(ymin, y);
            xmax = max(xmax, x);
            ymax = max(ymax, y);
        }

        void extendx(double r, double sintheta, double sign)
        {
            double x = r * sintheta * sign;
            xmin = min(xmin, x);
            xmax = max(xmax, x);
        }

        void extendy(double r, double sintheta, double sign)
        {
            double y = r * sintheta * sign;
            ymin = min(ymin, y);
            ymax = max(ymax, y);
        }
    };

    // initialize the equatorial bounds from the projected coordinates of the eight corner points
    EqBox eqbox(_rmin, _sinthetamin, _cosphimin, _sinphimin);
    eqbox.extend(_rmin, _sinthetamin, _cosphimax, _sinphimax);
    eqbox.extend(_rmin, _sinthetamax, _cosphimin, _sinphimin);
    eqbox.extend(_rmin, _sinthetamax, _cosphimax, _sinphimax);
    eqbox.extend(_rmax, _sinthetamin, _cosphimin, _sinphimin);
    eqbox.extend(_rmax, _sinthetamin, _cosphimax, _sinphimax);
    eqbox.extend(_rmax, _sinthetamax, _cosphimin, _sinphimin);
    eqbox.extend(_rmax, _sinthetamax, _cosphimax, _sinphimax);

    // extend with radial points if theta crosses the equatorial plane
    double thetacross = (_thetamin <= M_PI_2 && _thetamax >= M_PI_2);
    if (thetacross)
    {
        eqbox.extend(_rmax, 1., _cosphimin, _sinphimin);
        eqbox.extend(_rmax, 1., _cosphimax, _sinphimax);
    }

    // extend with radial points if phi crosses an axis (cannot cross negative x axis)
    double sintheta = thetacross ? 1. : max(_sinthetamin, _sinthetamax);
    if (_phimin <= -M_PI_2 && _phimax >= -M_PI_2) eqbox.extendy(_rmax, sintheta, -1.);  // negative y
    if (_phimin <= M_PI_2 && _phimax >= M_PI_2) eqbox.extendy(_rmax, sintheta, 1.);     // positive y
    if (_phimin <= 0. && _phimax >= 0.) eqbox.extendx(_rmax, sintheta, 1.);             // positive x

    return Box(eqbox.xmin, eqbox.ymin, zmin, eqbox.xmax, eqbox.ymax, zmax);
}

//////////////////////////////////////////////////////////////////////

double SphericalCell::intersection(Vec r, const Vec k) const
{
    // small value used to check for parallel directions
    constexpr double eps = 1e-12;

    // allocate room for the 10 possible intersections (2 per sphere, 2 per cone, and 1 per plane)
    // plus the starting position (which starts the first segment).
    // initialize the array to zeroes:
    //  - leave the first value at zero to represent the starting position
    //  - overwrite the other values for each intersection, or leave at zero if there is none
    enum { START, RMIN1, RMIN2, RMAX1, RMAX2, THETAMIN1, THETAMIN2, THETAMAX1, THETAMAX2, PHIMIN, PHIMAX, LEN };
    std::array<double, LEN> sv = {};

    // precalculate some properties of the ray
    double rk = Vec::dot(r, k);
    double r2 = r.norm2();
    double rz2 = r.z() * r.z();
    double kz2 = k.z() * k.z();
    double rkz = r.z() * k.z();

    // intersections with the rmin and rmax boundary spheres
    Quadratic::distinctSolutions(rk, r2 - _rmin * _rmin, sv[RMIN1], sv[RMIN2]);
    Quadratic::distinctSolutions(rk, r2 - _rmax * _rmax, sv[RMAX1], sv[RMAX2]);

    // intersections with the thetamin boundary cone
    {
        double cos2 = _costhetamin * _costhetamin;
        if (abs(cos2) >= eps)
        {
            double a = cos2 - kz2;
            double b = cos2 * rk - rkz;
            double c = cos2 * r2 - rz2;
            if (abs(a) >= eps)
            {
                // general case
                Quadratic::distinctSolutions(b / a, c / a, sv[THETAMIN1], sv[THETAMIN2]);
            }
            else
            {
                // ray parallel to cone
                if (abs(b) >= eps) sv[THETAMIN1] = -0.5 * c / b;
            }
        }
        else
        {
            // degenerate cone identical to xy-plane
            if (abs(k.z()) >= eps) sv[THETAMIN1] = -r.z() / k.z();
        }
    }

    // intersections with the thetamax boundary cone
    {
        double cos2 = _costhetamax * _costhetamax;
        if (abs(cos2) >= eps)
        {
            double a = cos2 - kz2;
            double b = cos2 * rk - rkz;
            double c = cos2 * r2 - rz2;
            if (abs(a) >= eps)
            {
                // general case
                Quadratic::distinctSolutions(b / a, c / a, sv[THETAMAX1], sv[THETAMAX2]);
            }
            else
            {
                // ray parallel to cone
                if (abs(b) >= eps) sv[THETAMAX1] = -0.5 * c / b;
            }
        }
        else
        {
            // degenerate cone identical to xy-plane
            if (abs(k.z()) >= eps) sv[THETAMAX1] = -r.z() / k.z();
        }
    }

    // intersection with meridional phimin plane
    {
        double q = k.x() * _sinphimin - k.y() * _cosphimin;
        if (abs(q) >= eps)
        {
            sv[PHIMIN] = -(r.x() * _sinphimin - r.y() * _cosphimin) / q;
        }
    }

    // intersection with meridional phimax plane
    {
        double q = k.x() * _sinphimax - k.y() * _cosphimax;
        if (abs(q) >= eps)
        {
            sv[PHIMAX] = -(r.x() * _sinphimax - r.y() * _cosphimax) / q;
        }
    }

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
