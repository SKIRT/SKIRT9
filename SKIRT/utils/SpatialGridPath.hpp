/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef SPATIALGRIDPATH_HPP
#define SPATIALGRIDPATH_HPP

#include "Direction.hpp"
#include "Position.hpp"
class Box;

//////////////////////////////////////////////////////////////////////

/** A SpatialGridPath object contains the geometric details of a path through a spatial grid. Given
    a spatial grid, i.e. some partition of space into cells, a starting position \f${\bf{r}}\f$ and
    a propagation direction \f${\bf{k}}\f$, one can calculate the path through the grid. A
    SpatialGridPath object maintains a record for each cell crossed by the path, called a \em
    segment. A segment stores the spatial cell index (so the cell can be identified in the grid),
    together with the physical path length \f$\Delta s\f$ covered within the cell, and the path
    length \f$s\f$ covered along the entire path up to the end of the cell.

    Updating the initial position and/or the direction of the path invalidates all segments in the
    path, but the segments are not automatically cleared. One should call the clear() function or
    the moveInside() function to do so.

    In addition to the geometric data, a SpatialGridPath object allows client code to store optical
    depth information with each path segment, as well as information on the interaction point along
    the path corresponding to a given optical depth. Specifically, each path segment can store
    either an extinction optical depth (for forced-scattering photon life cycles \em without
    explicit absorption) or both a scattering and absorption optical depth (for forced-scattering
    photon life cycles \em with explicit absorption). The interaction point information includes
    the spatial cell index, the cumulative distance, and the cumulative absorption optical depth
    (for the second type of photon cycle). */
class SpatialGridPath
{
public:
    // ------- Constructors; handling position and directions -------

    /** This constructor creates an empty path with the specified initial position and propagation
        direction. */
    SpatialGridPath(const Position& bfr, const Direction& bfk);

    /** This constructor creates an empty path with the initial position and propagation direction
        initialized to null values. After using this constructor, invoke the setPosition() and
        setDirection() functions to set these properties to appropriate values. */
    SpatialGridPath();

    /** This function sets the initial position of the path to a new value. */
    void setPosition(const Position& bfr) { _bfr = bfr; }

    /** This function sets the propagation direction along the path to a new value. */
    void setDirection(const Direction& bfk) { _bfk = bfk; }

    /** This function propagates the initial position of the path over a distance \f$s\f$. In other
        words, it updates the position from \f${\bf{r}}\f$ to \f${\bf{r}}+s\,{\bf{k}}\f$. */
    void propagatePosition(double s) { _bfr += s * _bfk; }

    /** This function returns the initial position of the path. */
    Position position() const { return _bfr; }

    /** This function returns the propagation direction along the path. */
    Direction direction() const { return _bfk; }

    // ------- Adding path segments -------

    /** This function removes all path segments, resulting in an empty path with the original
        initial position and propagation direction. */
    void clear();

    /** This function adds a segment in cell \f$m\f$ with length \f$\Delta s\f$ to the path.
        If \f$\Delta s\le 0\f$, the function does nothing. */
    void addSegment(int m, double ds);

    /** This function clears the path, adds any segments needed to move the initial position along
        the propagation direction (both specified in the constructor) inside a given box, and
        finally returns the resulting position. The small value specified by \em eps is added to
        the path length beyond the intersection point so that the final position is well inside the
        box, guarding against rounding errors. If the initial position is already inside the box,
        no segments are added. If the half-ray formed by the initial position and the propagation
        direction does not intersect the box, the function returns some arbitrary position outside
        the box. */
    Position moveInside(const Box& box, double eps);

    // ------- Working with path segments -------

    /** This data structure holds information about a given segment in the path. In addition to the
        spatial information (cell index, incremental distance, cumulative distance) managed by the
        SpatialGridPath class itself, segments can also store optical depth information managed by
        client code for forced-scattering photon life cycles. For photon life cycles \em without
        explicit absorption, just the extinction optical depth is stored. For life cycles \em with
        explicit absorption, the scattering and absorption optical depth are stored separately. */
    class Segment
    {
        int _m;               // cell index
        double _ds;           // distance within the cell
        double _s;            // cumulative distance until cell exit
        double _tauExtOrSca;  // cumulative extinction or scattering optical depth at cell exit
        double _tauAbs;       // cumulative absorption optical depth at cell exit or zero
    public:
        Segment(int m, double ds, double s) : _m(m), _ds(ds), _s(s), _tauExtOrSca{0.}, _tauAbs{0.} {}
        int m() const { return _m; }
        double ds() const { return _ds; }
        double s() const { return _s; }
        double tauExtOrSca() const { return _tauExtOrSca; }
        double tauAbs() const { return _tauAbs; }
        double tauExt() const { return _tauExtOrSca + _tauAbs; }
        void setOpticalDepth(double tauExt) { _tauExtOrSca = tauExt; }
        void setOpticalDepth(double tauSca, double tauAbs)
        {
            _tauExtOrSca = tauSca;
            _tauAbs = tauAbs;
        }
    };

    /** This function returns a read-only reference to the list of the current segments in the
        path. */
    const vector<Segment>& segments() const { return _segments; }

    /** This function returns a writable reference to the list of the current segments in the path.
        Because the segments are writable, the client code can update the stored optical depth
        information. */
    vector<Segment>& segments() { return _segments; }

    // ------- Handling the interaction point -------

    /** This function returns the optical depth corresponding to the end of the last path segment
        in the path, or zero if the path has no segments. The procedure uses the extinction or
        scattering optical depth depending on which one has been stored in the segments by the
        client code.

        The function assumes that both the geometric and optical depth information for the path
        have been set; if this is not the case, the behavior is undefined. */
    double totalOpticalDepth() const;

    /** This function determines the interaction point along the path corresponding to the
        specified interaction optical depth, and stores relevant information about it in data
        members for later retrieval through the interactionCellIndex(), interactionDistance() and
        interactionOpticalDepth() functions.

        Specifically, the function first determines the path segment for which the exit optical
        depth becomes larger than the specified interaction optical depth. The procedure uses the
        extinction or scattering optical depth depending on which one has been stored in the
        segments by the client code. The function then calculates the distance covered along the
        path until the interaction optical depth has been reached by linear interpolation within
        the interacting cell. If applicable, it similarly interpolates and the absorption optical
        depth at the interaction point. Using linear interpolation is equivalent to assuming
        exponential behavior of the extinction with distance within the cell. Finally, the function
        stores the index of the spatial cell corresponding to the interacting segment plus the
        interpolated values in data members for later retrieval.

        The function assumes that both the geometric and optical depth information for the path
        have been set; if this is not the case, the behavior is undefined. */
    void findInteractionPoint(double tauinteract);

    /** This function stores the specified spatial cell index, distance and cumulative absorption
        optical depth to the initial position of the interaction point for later retrieval through
        the interactionCellIndex(), interactionDistance(), and interactionOpticalDepth() functions.
        */
    void setInteractionPoint(int m, double s, double tauAbs = 0.);

    /** This function returns the spatial cell index corresponding to the interaction point most
        recently set by one of the findInteractionPoint() or setInteractionPoint() functions, or -1
        if these functions have never been called or if there was no interaction point within the
        path. */
    int interactionCellIndex() const { return _interactionCellIndex; }

    /** This function returns the distance along the path from its initial position to the
        interaction point most recently set by one of the findInteractionPoint() or
        setInteractionPoint() functions, or zero if these functions have never been called or if
        there was no interaction point within the path. */
    double interactionDistance() const { return _interactionDistance; }

    /** This function returns the cumulative absorption optical depth along the path from its
        initial position to the interaction point most recently set by one of the
        findInteractionPoint() or setInteractionPoint() functions, or zero if these functions has
        never been called. */
    double interactionOpticalDepth() const { return _interactionOpticalDepth; }

    // ------- Data members -------
private:
    Position _bfr;
    Direction _bfk;
    vector<Segment> _segments;
    double _s{0.};
    int _interactionCellIndex{-1};
    double _interactionDistance{0.};
    double _interactionOpticalDepth{0.};
};

//////////////////////////////////////////////////////////////////////

#endif
