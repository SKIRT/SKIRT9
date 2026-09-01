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
