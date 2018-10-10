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

double SpatialGridPath::totalOpticalDepth()
{
    return !_segments.empty() ? _segments.back().tau : 0.;
}

////////////////////////////////////////////////////////////////////

void SpatialGridPath::findInteractionPoint(double tau)
{
    // we can't handle an empty path
    if (_segments.empty())
    {
        _interactionCellIndex = -1;
        _interactionDistance = 0.;
    }
    else
    {
        // find a pointer to the first segment that has an exit optical depth larger than or equal to the given value,
        // or a pointer beyond the list if no such element is found
        auto seg = std::lower_bound(_segments.cbegin(), _segments.cend(), tau,
                                    [](const Segment& seg, double t) {return seg.tau < t; });

        // if we find the first segment, interpolate with the path's entry point
        if (seg == _segments.cbegin())
        {
            _interactionCellIndex = seg->m;
            _interactionDistance = NR::interpolateLinLin(tau, 0., seg->tau, 0., seg->s);
        }

        // if we find some other segment, interpolate with the previous segment
        else if (seg < _segments.cend())
        {
            _interactionCellIndex = seg->m;
            _interactionDistance = NR::interpolateLinLin(tau, (seg-1)->tau, seg->tau, (seg-1)->s, seg->s);
        }

        // if we are beyond the last segment, just use the last segment (i.e. assume this is a numerical inaccuracy)
        else
        {
            _interactionCellIndex = (seg-1)->m;
            _interactionDistance = (seg-1)->s;
        }
    }
}

////////////////////////////////////////////////////////////////////
