/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "SpatialGridPath.hpp"
#include "Box.hpp"
#include "NR.hpp"

////////////////////////////////////////////////////////////////////

namespace
{
    const int INITIAL_CAPACITY = 1000;
}

////////////////////////////////////////////////////////////////////

SpatialGridPath::SpatialGridPath(const Position& bfr, const Direction& bfk)
    : _bfr(bfr), _bfk(bfk)
{
    _segments.reserve(INITIAL_CAPACITY);
}

////////////////////////////////////////////////////////////////////

SpatialGridPath::SpatialGridPath()
{
    _segments.reserve(INITIAL_CAPACITY);
}

////////////////////////////////////////////////////////////////////

void SpatialGridPath::clear()
{
    _segments.clear();
}

////////////////////////////////////////////////////////////////////

void SpatialGridPath::addSegment(int m, double ds)
{
    if (ds>0)
    {
        double s = !_segments.empty() ? _segments.back().s : 0.;
        _segments.push_back(Segment{m, ds, s+ds, 0.});
    }
}

////////////////////////////////////////////////////////////////////

Position SpatialGridPath::moveInside(const Box& box, double eps)
{
    // a position that is certainly not inside any box
    static const Position OUTSIDE(std::numeric_limits<double>::infinity(),
                                  std::numeric_limits<double>::infinity(),
                                  std::numeric_limits<double>::infinity());

    // clear the path
    clear();

    // initial position and direction
    double kx,ky,kz;
    _bfk.cartesian(kx,ky,kz);
    double rx,ry,rz;
    _bfr.cartesian(rx,ry,rz);

    // --> x direction
    if (rx <= box.xmin())
    {
        if (kx <= 0.0) return OUTSIDE;
        else
        {
            double ds = (box.xmin()-rx)/kx;
            addSegment(-1,ds);
            rx = box.xmin() + eps;
            ry += ky*ds;
            rz += kz*ds;
        }
    }
    else if (rx >= box.xmax())
    {
        if (kx >= 0.0) return OUTSIDE;
        else
        {
            double ds = (box.xmax()-rx)/kx;
            addSegment(-1,ds);
            rx = box.xmax() - eps;
            ry += ky*ds;
            rz += kz*ds;
        }
    }

    // --> y direction
    if (ry <= box.ymin())
    {
        if (ky <= 0.0) return OUTSIDE;
        else
        {
            double ds = (box.ymin()-ry)/ky;
            addSegment(-1,ds);
            rx += kx*ds;
            ry = box.ymin() + eps;
            rz += kz*ds;
        }
    }
    else if (ry >= box.ymax())
    {
        if (ky >= 0.0) return OUTSIDE;
        else
        {
            double ds = (box.ymax()-ry)/ky;
            addSegment(-1,ds);
            rx += kx*ds;
            ry = box.ymax() - eps;
            rz += kz*ds;
        }
    }

    // --> z direction
    if (rz <= box.zmin())
    {
        if (kz <= 0.0) return OUTSIDE;
        else
        {
            double ds = (box.zmin()-rz)/kz;
            addSegment(-1,ds);
            rx += kx*ds;
            ry += ky*ds;
            rz = box.zmin() + eps;
        }
    }
    else if (rz >= box.zmax())
    {
        if (kz >= 0.0) return OUTSIDE;
        else
        {
            double ds = (box.zmax()-rz)/kz;
            addSegment(-1,ds);
            rx += kx*ds;
            ry += ky*ds;
            rz = box.zmax() - eps;
        }
    }

    // the position should now be just inside the box; although in rare cases, it may be still be outside!
    return Position(rx,ry,rz);
}

////////////////////////////////////////////////////////////////////

double SpatialGridPath::escapeExtinctionFactor()
{
    return !_segments.empty() ? _segments.back().zeta : 1.;
}

////////////////////////////////////////////////////////////////////

void SpatialGridPath::findInteractionPoint(double zeta)
{
    _interactionCellIndex = -1;
    _interactionDistance = 0.;

    // we can't handle an empty path, nor a specified extinction factor of one (or more)
    if (!_segments.empty() && zeta<1)
    {
        // find a pointer to the first segment that has an extinction factor smaller than or equal to the given value,
        // or a pointer beyond the list if no such element is found
        auto p = std::upper_bound(_segments.cbegin(), _segments.cend(), zeta,
                                  [](double z, const Segment& seg) {return z >= seg.zeta; });

        // if we find the first segment, interpolate with the extinction factor for the path's entry point (i.e. 1.)
        if (p == _segments.cbegin())
        {
            _interactionCellIndex = _segments[0].m;
            _interactionDistance = NR::interpolateLogLin(zeta, 1., _segments[0].zeta, 0., _segments[0].s);
        }

        // if we find some other segment, interpolate with the previous segment
        else if (p < _segments.cend())
        {
            auto i = p - _segments.cbegin();
            _interactionCellIndex = _segments[i].m;
            _interactionDistance = NR::interpolateLogLin(zeta, _segments[i-1].zeta, _segments[i].zeta,
                                                               _segments[i-1].s, _segments[i].s);
        }

        // if we are beyond the last segment, just use the last segment (i.e. assume this is a numerical inaccuracy)
        else
        {
            _interactionCellIndex = _segments.back().m;
            _interactionDistance = _segments.back().s;
        }
    }
}

////////////////////////////////////////////////////////////////////
