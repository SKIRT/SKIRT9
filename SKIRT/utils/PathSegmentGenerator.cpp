/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observato_ry(), Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "PathSegmentGenerator.hpp"
#include "Box.hpp"

////////////////////////////////////////////////////////////////////

bool PathSegmentGenerator::moveInside(const Box& box, double eps)
{
    // initialize to empty segment with zero length
    setEmptySegment();
    _state = State::Outside;

    // keep track of cumulative length of all subsegments
    double cumds = 0.;

    // --> x direction
    if (_rx <= box.xmin())
    {
        if (_kx <= 0.0)
            return false;
        else
        {
            double ds = (box.xmin() - _rx) / _kx;
            _rx = box.xmin() + eps;
            _ry += _ky * ds;
            _rz += _kz * ds;
            cumds += ds;
        }
    }
    else if (_rx >= box.xmax())
    {
        if (_kx >= 0.0)
            return false;
        else
        {
            double ds = (box.xmax() - _rx) / _kx;
            _rx = box.xmax() - eps;
            _ry += _ky * ds;
            _rz += _kz * ds;
            cumds += ds;
        }
    }

    // --> y direction
    if (_ry <= box.ymin())
    {
        if (_ky <= 0.0)
            return false;
        else
        {
            double ds = (box.ymin() - _ry) / _ky;
            _rx += _kx * ds;
            _ry = box.ymin() + eps;
            _rz += _kz * ds;
            cumds += ds;
        }
    }
    else if (_ry >= box.ymax())
    {
        if (_ky >= 0.0)
            return false;
        else
        {
            double ds = (box.ymax() - _ry) / _ky;
            _rx += _kx * ds;
            _ry = box.ymax() - eps;
            _rz += _kz * ds;
            cumds += ds;
        }
    }

    // --> z direction
    if (_rz <= box.zmin())
    {
        if (_kz <= 0.0)
            return false;
        else
        {
            double ds = (box.zmin() - _rz) / _kz;
            _rx += _kx * ds;
            _ry += _ky * ds;
            _rz = box.zmin() + eps;
            cumds += ds;
        }
    }
    else if (_rz >= box.zmax())
    {
        if (_kz >= 0.0)
            return false;
        else
        {
            double ds = (box.zmax() - _rz) / _kz;
            _rx += _kx * ds;
            _ry += _ky * ds;
            _rz = box.zmax() - eps;
            cumds += ds;
        }
    }

    // in rare border cases, the position can still be outside of the box
    if (!box.contains(r())) return false;

    // return the empty segment with the cumulative length
    // (which might be zero if the position was inside the box to begin with)
    setEmptySegment(cumds);
    _state = State::Inside;
    return true;
}

////////////////////////////////////////////////////////////////////
