/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef PATHSEGMENTGENERATOR_HPP
#define PATHSEGMENTGENERATOR_HPP

#include "SpatialGridPath.hpp"

//////////////////////////////////////////////////////////////////////

/** PathSegmentGenerator is the abstract base class for classes that calculate and return the
    segments of a path through a spatial grid segment by segment. This base class offers the public
    interface for obtaining the segments one by one after initializing the path with a given
    starting position and direction. A subclass must be provided for each type of spatial grid.
    Usually, this subclass is implemented as a private class inside the corresponding SpatialGrid
    subclass.

    In addition to the public interface, the PathSegmentGenerator class also offers facilities for
    use by subclasses. These include functions to access the position and direction of the path, to
    set the cell index and path length for the next segment, and to help track the state of the
    generator. */
class PathSegmentGenerator
{
    // ------- Constructing and destructing -------

public:
    /** The constructor creates a path segment generator initialized to an invalid path. To start
        path segment generation, call the start() function. */
    PathSegmentGenerator() {}

    /** The destructor does what it is expected to do. */
    virtual ~PathSegmentGenerator() {}

    // ------- Generating and retrieving path segments -------

public:
    /** This function initializes path segment generation for the starting position \f${\bf{r}}\f$
        and the direction \f${\bf{k}}\f$ specified by the SpatialGridPath instance passed as an
        argument. The function \em must be called before calling the next() function for the first
        time, or to re-initialize the generator for a fresh path. */
    void start(const SpatialGridPath* path)
    {
        _state = State::Unknown;
        path->position().cartesian(_rx, _ry, _rz);
        path->direction().cartesian(_kx, _ky, _kz);
    }

    /** This function initializes path segment generation for the specified starting position
        \f${\bf{r}}\f$ and direction \f${\bf{k}}\f$. The function \em must be called before calling
        the next() function for the first time, or to re-initialize the generator for a fresh path.
        */
    void start(Position bfr, Direction bfk)
    {
        _state = State::Unknown;
        bfr.cartesian(_rx, _ry, _rz);
        bfk.cartesian(_kx, _ky, _kz);
    }

    /** This function calculates the next path segment and stores its cell index and path length in
        data members that can be accessed through the m() and ds() functions. It should be called
        only after the path has been intialized through the start() function. The next() function
        returns true if a path segment is available and false if there are no more segments in the
        path. In the latter case, the values returned by the m() and ds() functions are undefined.

        This function must be implemented by each subclass. */
    virtual bool next() = 0;

    /** This function returns the index of the spatial cell crossed by the current path segment.
        The value is defined only after a call to the next() function returned true. */
    int m() { return _m; }

    /** This function returns the length of the current path segment inside the spatial cell
        crossed by the segment. The value is defined only after a call to the next() function
        returned true. */
    double ds() { return _ds; }

    // ------- Accessing internal state - for use by subclasses -------

protected:
    /** This enumeration lists the states that are passed through by most path segment generators.
        The start() function sets the state to Unknown. */
    enum class State { Unknown, Inside, Outside };

    /** This function returns the current state of the generator. */
    State state() const { return _state; }

    /** This function sets the current state of the generator to the specified value. */
    void setState(State state) { _state = state; }

    /** This function returns the path's position. This value is initialized by the start()
        function and can be updated from a subclass through various functions offered by this base
        class. */
    Position r() { return Position(_rx, _ry, _rz); }

    /** This function returns the x component of the path's position. This value is initialized by
        the start() function and can be updated from a subclass through various functions offered
        by this base class. */
    double rx() const { return _rx; }

    /** This function returns the y component of the path's position. This value is initialized by
        the start() function and can be updated from a subclass through various functions offered
        by this base class. */
    double ry() const { return _ry; }

    /** This function returns the y component of the path's position. This value is initialized by
        the start() function and can be updated from a subclass through various functions offered
        by this base class. */
    double rz() const { return _rz; }

    /** This function returns the path's direction. This value is initialized by the start()
        function and remains immutable during segment generation for the path. */
    Direction k() { return Direction(_kx, _ky, _kz); }

    /** This function returns the x component of the path's direction. This value is initialized by
        the start() function and remains immutable during segment generation for the path. */
    double kx() const { return _kx; }

    /** This function returns the y component of the path's direction. This value is initialized by
        the start() function and remains immutable during segment generation for the path. */
    double ky() const { return _ky; }

    /** This function returns the z component of the path's direction. This value is initialized by
        the start() function and remains immutable during segment generation for the path. */
    double kz() const { return _kz; }

    /** This function sets the x component of the path's position to the specified value. */
    void setrx(double rx) { _rx = rx; }

    /** This function sets the y component of the path's position to the specified value. */
    void setry(double ry) { _ry = ry; }

    /** This function sets the z component of the path's position to the specified value. */
    void setrz(double rz) { _rz = rz; }

    /** This function advances the path's position along the path's direction for the specified
        distance. */
    void propagater(double ds)
    {
        _rx += _kx * ds;
        _ry += _ky * ds;
        _rz += _kz * ds;
    }

    /** This function advances each component of the path's position to the next representable
        double-precision flaoting point number in the path's direction. */
    void propagateToNextAfter()
    {
        _rx = std::nextafter(_rx, (_kx < 0.) ? -DBL_MAX : DBL_MAX);
        _ry = std::nextafter(_ry, (_ky < 0.) ? -DBL_MAX : DBL_MAX);
        _rz = std::nextafter(_rz, (_kz < 0.) ? -DBL_MAX : DBL_MAX);
    }

    /** This function advances the x component of the path's position along the path's direction
        for the specified distance. */
    void propagaterx(double ds) { _rx += _kx * ds; }

    /** This function advances the y component of the path's position along the path's direction
        for the specified distance. */
    void propagatery(double ds) { _ry += _ky * ds; }

    /** This function advances the z component of the path's position along the path's direction
        for the specified distance. */
    void propagaterz(double ds) { _rz += _kz * ds; }

    /** This function sets the information for the current path segment to a cell number of -1,
        indicating a path segment outside of the grid, and the specified segment length. If the
        argument is omitted, a length of zero is assumed. */
    void setEmptySegment(double ds = 0.)
    {
        _m = -1;
        _ds = ds;
    }

    /** This function sets the information for the current path segment to the specified cell
        number and segment length. */
    void setSegment(int m, double ds)
    {
        _m = m;
        _ds = ds;
    }

    /** This function calculates the segment needed to advance the path along its direction from
        its current position to the closest position inside the specified box (a cuboid aligned
        with the coordinate axes).

        If the path does not intersect the box, the function returns false. In that case, the state
        is set to Outside and no segment is available (the path position and segment values are
        undefined). If the original position is inside the box or the path does intersect the box,
        the function returns true. In that case, the state is set to Inside, the position is
        adjusted to be inside the box if necessary, the segment cell index is set to -1 indicating
        a segment outside of the box, and the segment length is set to the distance between the new
        and original positions. If the original position was inside the box, the segment length is
        set to zero. */
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
