/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef SPATIALGRIDSEGMENTGENERATOR_HPP
#define SPATIALGRIDSEGMENTGENERATOR_HPP

#include "SpatialGridPath.hpp"

//////////////////////////////////////////////////////////////////////

/** A SpatialGridSegmentGenerator object -- TO DO. */
class SpatialGridSegmentGenerator
{
    // ------- Constructing and destructing -------

public:
    /** TO DO. */
    SpatialGridSegmentGenerator(const SpatialGridPath* path)
        : _rx{path->position().x()}, _ry{path->position().y()}, _rz{path->position().z()}, _kx{path->direction().x()},
          _ky{path->direction().y()}, _kz{path->direction().z()}
    {}

    /** TO DO. */
    virtual ~SpatialGridSegmentGenerator() {}

    // ------- Generating and retrieving path segments -------

public:
    /** True if segment is available; false if no more segments. */
    virtual bool next() = 0;

    /** TO DO. */
    int m() { return _m; }

    /** TO DO. */
    double ds() { return _ds; }

    // ------- Accessing internal state - for use by subclasses -------

protected:
    Position r() { return Position(_rx, _ry, _rz); }

    double rx() { return _rx; }
    double ry() { return _ry; }
    double rz() { return _rz; }

    double kx() { return _kx; }
    double ky() { return _ky; }
    double kz() { return _kz; }

    void setrx(double rx) { _rx = rx; }
    void setry(double ry) { _ry = ry; }
    void setrz(double rz) { _rz = rz; }

    void setoutsidepos() { _rx = _ry = _rz = std::numeric_limits<double>::infinity(); }

    void propagatex() { _rx += _kx * _ds; }
    void propagatey() { _ry += _ky * _ds; }
    void propagatez() { _rz += _kz * _ds; }

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

    /** sets the segment and adjusts the position */
    void moveInside(const Box& box, double eps);

    // ------- Data members -------

private:
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
