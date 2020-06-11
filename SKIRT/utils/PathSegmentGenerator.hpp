/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef PATHSEGMENTGENERATOR_HPP
#define PATHSEGMENTGENERATOR_HPP

#include "SpatialGridPath.hpp"

//////////////////////////////////////////////////////////////////////

/** A PathSegmentGenerator object -- TO DO. */
class PathSegmentGenerator
{
    // ------- Constructing and destructing -------

public:
    /** TO DO. */
    PathSegmentGenerator() {}

    /** TO DO. */
    virtual ~PathSegmentGenerator() {}

    // ------- Generating and retrieving path segments -------

public:
    /** TO DO. */
    void start(const SpatialGridPath* path)
    {
        _state = State::Unknown;
        path->position().cartesian(_rx, _ry, _rz);
        path->direction().cartesian(_kx, _ky, _kz);
    }

    /** True if segment is available; false if no more segments. */
    virtual bool next() = 0;

    /** TO DO. */
    int m() { return _m; }

    /** TO DO. */
    double ds() { return _ds; }

    // ------- Accessing internal state - for use by subclasses -------

protected:
    enum class State { Unknown, Inside, Outside };
    State state() const { return _state; }
    void setState(State state) { _state = state; }

    Position r() { return Position(_rx, _ry, _rz); }
    Direction k() { return Direction(_kx, _ky, _kz); }

    /** TO DO: streamline these functions and their usage by the various grids) */
    double rx() const { return _rx; }
    double ry() const { return _ry; }
    double rz() const { return _rz; }

    double kx() const { return _kx; }
    double ky() const { return _ky; }
    double kz() const { return _kz; }

    void setrx(double rx) { _rx = rx; }
    void setry(double ry) { _ry = ry; }
    void setrz(double rz) { _rz = rz; }

    void propagater(double ds)
    {
        _rx += _kx * ds;
        _ry += _ky * ds;
        _rz += _kz * ds;
    }
    void propagaterx(double ds) { _rx += _kx * ds; }
    void propagatery(double ds) { _ry += _ky * ds; }
    void propagaterz(double ds) { _rz += _kz * ds; }

    void setSegment(int m, double ds)
    {
        _m = m;
        _ds = ds;
    }

    void setEmptySegment(double ds = 0.)
    {
        _m = -1;
        _ds = ds;
    }

    /** sets the segment and adjusts the position and state; returns true if position is now inside */
    bool moveInside(const Box& box, double eps);

    // ------- Data members -------

private:
    State _state{State::Unknown};
    double _rx{0.};
    double _ry{0.};
    double _rz{0.};
    double _kx{0.};
    double _ky{0.};
    double _kz{0.};
    double _ds{0.};
    int _m{-1};
};

//////////////////////////////////////////////////////////////////////

#endif
