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

SpatialGridPath::SpatialGridPath(const Position& bfr, const Direction& bfk) : _bfr(bfr), _bfk(bfk)
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
    _s = 0.;
}

////////////////////////////////////////////////////////////////////

void SpatialGridPath::addSegment(int m, double ds)
{
    if (ds > 0.)
    {
        _s += ds;
        _segments.emplace_back(m, ds, _s);
    }
}

////////////////////////////////////////////////////////////////////

Position SpatialGridPath::moveInside(const Box& box, double eps)
{
    // a position that is certainly not inside any box
    static const Position OUTSIDE(std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity(),
                                  std::numeric_limits<double>::infinity());

    // clear the path
    clear();

    // initial position and direction
    double kx, ky, kz;
    _bfk.cartesian(kx, ky, kz);
    double rx, ry, rz;
    _bfr.cartesian(rx, ry, rz);

    // --> x direction
    if (rx <= box.xmin())
    {
        if (kx <= 0.0)
            return OUTSIDE;
        else
        {
            double ds = (box.xmin() - rx) / kx;
            addSegment(-1, ds);
            rx = box.xmin() + eps;
            ry += ky * ds;
            rz += kz * ds;
        }
    }
    else if (rx >= box.xmax())
    {
        if (kx >= 0.0)
            return OUTSIDE;
        else
        {
            double ds = (box.xmax() - rx) / kx;
            addSegment(-1, ds);
            rx = box.xmax() - eps;
            ry += ky * ds;
            rz += kz * ds;
        }
    }

    // --> y direction
    if (ry <= box.ymin())
    {
        if (ky <= 0.0)
            return OUTSIDE;
        else
        {
            double ds = (box.ymin() - ry) / ky;
            addSegment(-1, ds);
            rx += kx * ds;
            ry = box.ymin() + eps;
            rz += kz * ds;
        }
    }
    else if (ry >= box.ymax())
    {
        if (ky >= 0.0)
            return OUTSIDE;
        else
        {
            double ds = (box.ymax() - ry) / ky;
            addSegment(-1, ds);
            rx += kx * ds;
            ry = box.ymax() - eps;
            rz += kz * ds;
        }
    }

    // --> z direction
    if (rz <= box.zmin())
    {
        if (kz <= 0.0)
            return OUTSIDE;
        else
        {
            double ds = (box.zmin() - rz) / kz;
            addSegment(-1, ds);
            rx += kx * ds;
            ry += ky * ds;
            rz = box.zmin() + eps;
        }
    }
    else if (rz >= box.zmax())
    {
        if (kz >= 0.0)
            return OUTSIDE;
        else
        {
            double ds = (box.zmax() - rz) / kz;
            addSegment(-1, ds);
            rx += kx * ds;
            ry += ky * ds;
            rz = box.zmax() - eps;
        }
    }

    // the position should now be just inside the box; although in rare cases, it may be still be outside!
    return Position(rx, ry, rz);
}

////////////////////////////////////////////////////////////////////

double SpatialGridPath::totalOpticalDepth() const
{
    return !_segments.empty() ? _segments.back().tauExtOrSca() : 0.;
}

////////////////////////////////////////////////////////////////////

void SpatialGridPath::findInteractionPoint(double tauinteract)
{
    // we can't handle an empty path
    if (_segments.empty())
    {
        _interactionCellIndex = -1;
        _interactionDistance = 0.;
        _interactionOpticalDepth = 0.;
    }
    else
    {
        // find a pointer to the first segment that has an exit optical depth strictly larger than the given value,
        // (so that we never select an empty segment) or a pointer beyond the list if no such element is found
        auto seg = std::upper_bound(_segments.cbegin(), _segments.cend(), tauinteract,
                                    [](double t, const Segment& seg) { return t < seg.tauExtOrSca(); });

        // if we find the first segment, interpolate with the path's entry point
        if (seg == _segments.cbegin())
        {
            _interactionCellIndex = seg->m();
            _interactionDistance = NR::interpolateLinLin(tauinteract, 0., seg->tauExtOrSca(), 0., seg->s());
            _interactionOpticalDepth = NR::interpolateLinLin(tauinteract, 0., seg->tauExtOrSca(), 0., seg->tauAbs());
        }

        // if we find some other segment, interpolate with the previous segment
        else if (seg < _segments.cend())
        {
            _interactionCellIndex = seg->m();
            _interactionDistance = NR::interpolateLinLin(tauinteract, (seg - 1)->tauExtOrSca(), seg->tauExtOrSca(),
                                                         (seg - 1)->s(), seg->s());
            _interactionOpticalDepth = NR::interpolateLinLin(tauinteract, (seg - 1)->tauExtOrSca(), seg->tauExtOrSca(),
                                                             (seg - 1)->tauAbs(), seg->tauAbs());
        }

        // if we are precisely at or beyond the exit optical depth of the last segment, just use the last segment
        else
        {
            _interactionCellIndex = (seg - 1)->m();
            _interactionDistance = (seg - 1)->s();
            _interactionOpticalDepth = (seg - 1)->tauAbs();
        }
    }
}

////////////////////////////////////////////////////////////////////

void SpatialGridPath::setInteractionPoint(int m, double s, double tauAbs)
{
    _interactionCellIndex = m;
    _interactionDistance = s;
    _interactionOpticalDepth = tauAbs;
}

////////////////////////////////////////////////////////////////////
