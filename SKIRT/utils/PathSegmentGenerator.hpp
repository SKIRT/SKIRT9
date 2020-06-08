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
    PathSegmentGenerator(const SpatialGridPath* path)
        : _rx{path->position().x()}, _ry{path->position().y()}, _rz{path->position().z()}, _kx{path->direction().x()},
          _ky{path->direction().y()}, _kz{path->direction().z()}
    {}

    /** TO DO. */
    virtual ~PathSegmentGenerator() {}

    // ------- Generating and retrieving path segments -------

public:
    /** True if segment is available; false if no more segments. */
    virtual bool next() = 0;

    /** TO DO. */
    int m() const { return _m; }

    /** TO DO. */
    double ds() const { return _ds; }

    // ------- Accessing internal state - for use by subclasses -------

protected:
    Position r() { return Position(_rx, _ry, _rz); }

    double rx() const { return _rx; }
    double ry() const { return _ry; }
    double rz() const { return _rz; }

    double kx() const { return _kx; }
    double ky() const { return _ky; }
    double kz() const { return _kz; }

    void setrx(double rx) { _rx = rx; }
    void setry(double ry) { _ry = ry; }
    void setrz(double rz) { _rz = rz; }

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

    /** sets the segment and adjusts the position; returns true if position is now inside */
    /** TO DO. This function clears the path, adds any segments needed to move the initial position along
        the propagation direction (both specified in the constructor) inside a given box, and
        finally returns the resulting position. The small value specified by \em eps is added to
        the path length beyond the intersection point so that the final position is well inside the
        box, guarding against rounding errors. If the initial position is already inside the box,
        no segments are added. If the half-ray formed by the initial position and the propagation
        direction does not intersect the box, the function returns some arbitrary position outside
        the box. */
    bool moveInside(const Box& box, double eps);

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
