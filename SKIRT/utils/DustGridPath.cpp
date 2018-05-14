/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "DustGridPath.hpp"
#include "Box.hpp"
#include "NR.hpp"

//////////////////////////////////////////////////////////////////////

namespace
{
    const int INITIAL_CAPACITY = 1000;
}

//////////////////////////////////////////////////////////////////////

DustGridPath::DustGridPath(const Position& bfr, const Direction& bfk)
    : _bfr(bfr), _bfk(bfk), _s(0)
{
    _v.reserve(INITIAL_CAPACITY);
}

//////////////////////////////////////////////////////////////////////

DustGridPath::DustGridPath()
    : _s(0)
{
    _v.reserve(INITIAL_CAPACITY);
}

//////////////////////////////////////////////////////////////////////

void DustGridPath::clear()
{
    _s = 0;
    _v.clear();
}

//////////////////////////////////////////////////////////////////////

void DustGridPath::addSegment(int m, double ds)
{
    if (ds>0)
    {
        _s += ds;
        _v.push_back(Segment{m,ds,_s,0,0});
    }
}

//////////////////////////////////////////////////////////////////////

Position DustGridPath::moveInside(const Box& box, double eps)
{
    // a position that is certainly not inside any box
    static const Position OUTSIDE(std::numeric_limits<double>::infinity(),
                                  std::numeric_limits<double>::infinity(),
                                  std::numeric_limits<double>::infinity());

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

//////////////////////////////////////////////////////////////////////

double DustGridPath::tau() const
{
    int N = _v.size();
    return N ? _v[N-1].tau : 0;
}

//////////////////////////////////////////////////////////////////////

double DustGridPath::pathLength(double tau) const
{
    int N = _v.size();
    if (N>0 && tau>0)
    {
        int i = NR::locate(_v,Segment{0,0,0,0,tau});
        if (i<0) return NR::interpolateLinLin(tau, 0, _v[0].tau, 0, _v[0].s);
        if (i<N-1) return NR::interpolateLinLin(tau, _v[i].tau, _v[i+1].tau, _v[i].s, _v[i+1].s);
        return _v[N-1].s;
    }
    return 0;
}

//////////////////////////////////////////////////////////////////////
