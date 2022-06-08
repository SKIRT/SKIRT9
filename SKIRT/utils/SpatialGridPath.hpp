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

    In addition, a SpatialGridPath object allows storing a (cumulative) optical depth \f$tau\f$
    with each path segment. Before this non-geometric information is used, it must be set
    explicitly by a client class after the segments of the path have been calculated. As a
    convenience to client classes, the SpatialGridPath class offers some operations on this extra
    information, such as locating the interaction point along the path corresponding to a given
    optical depth.

    Updating the initial position and/or the direction of the path invalidates all segments in the
    path, but the segments are not automatically cleared. One should call the clear() function or
    the moveInside() function to do so. */
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

    // ------- Retrieving path segments -------

    // basic data structure holding information about a given segment in the path
    struct Segment
    {
        int m;          // cell index
        double ds;      // distance within the cell
        double s;       // cumulative distance until cell exit
        double tau;     // cumulative optical depth at cell exit
                        //    regular photon cycle: extinction (scattering + absorption) optical depth
                        //    with explicit absorption: scattering optical depth
        double tauAbs;  // cumulative absorption optical depth at cell exit (only with explicit absorption)
    };

    /** This function returns (a read-only reference to) a list of the current segments in the
        path. Each Segment object holds the following public data members: the index \em m of the
        spatial cell being represented, the distance \f$\Delta s\f$ covered by the path inside the
        cell, the cumulative distance \f$s\f$ covered by the path from its initial position to the
        exit point from the cell, and the cumulative cumulative optical depth \f$\tau\f$ at the
        exit point (or zero if this information has not been set for the path). */
    const vector<Segment>& segments() const { return _segments; }

    // ------- Handling optical depth -------

    /** This function sets the extinction (scattering plus absorption) optical depth corresponding
        to the end of the path segment with zero-based index \f$i\f$ to the specified value. The
        optical depth values for consecutive segments must form a non-descending sequence, i.e. \f[
        0\leq\tau_0 \leq\tau_1 \leq\,\dots\, \leq\tau_{N-1}\f] where \f$N\f$ is the number of
        segments in the path.

        This (non-geometric) information can be stored here as a convenience to client classes. If
        the index is out of range, undefined behavior results. */
    void setOpticalDepth(int i, double tau) { _segments[i].tau = tau; }

    /** This function sets the scattering and absorption optical
        depths corresponding to the end of the path segment with zero-based index \f$i\f$ to the
        specified values. This function is used in photon cycles with explicit absorption. The
        scattering optical depth values for consecutive segments must form a non-descending
        sequence, but this restriction does not hold for the absorption optical depth.

        This (non-geometric) information can be stored here as a convenience to client classes. If
        the index is out of range, undefined behavior results. */
    void setOpticalDepth(int i, double tauSca, double tauAbs)
    {
        _segments[i].tau = tauSca;
        _segments[i].tauAbs = tauAbs;
    }

    /** This function returns the optical depth corresponding to the end of the last path segment
        in the path, or zero if the path has no segments. The function assumes that both the
        geometric and optical depth information for the path have been set; if this is not the
        case, the behavior is undefined. */
    double totalOpticalDepth() const;

    /** This function determines the interaction point along the path corresponding to the
        specified optical depth, and stores relevant information about it in data members for later
        retrieval through the interactionCellIndex() and interactionDistance() functions.
        Specifically, the function first determines the path segment for which the exit optical
        depth becomes smaller than the specified optical depth. It then stores the index of the
        spatial cell corresponding to the interacting segment. Finally, it calculates and stores
        the distance covered along the path until the specified optical depth has been reached by
        linear interpolation within the interacting cell. Using linear interpolation is equivalent
        to assuming exponential behavior of the extinction with distance within the cell.

        The function assumes that both the geometric and optical depth information for the path
        have been set; if this is not the case, the behavior is undefined. */
    void findInteractionPoint(double tau);

    /** This function determines the interaction point along the path corresponding to the
        specified optical depth, and stores relevant information about it in data members for later
        retrieval through the interactionCellIndex() and interactionDistance() functions.
        Specifically, the function first determines the path segment for which the exit optical
        depth becomes smaller than the specified optical depth. It then stores the index of the
        spatial cell corresponding to the interacting segment. Finally, it calculates and stores
        the distance covered along the path until the specified optical depth has been reached by
        linear interpolation within the interacting cell. Using linear interpolation is equivalent
        to assuming exponential behavior of the extinction with distance within the cell.

        The function assumes that both the geometric and optical depth information for the path
        have been set; if this is not the case, the behavior is undefined. */
    void findExplicitAbsorptionInteractionPoint(double tau);

    /** This function stores the specified spatial cell index and distance to the initial position
        of the interaction point for later retrieval through the interactionCellIndex() and
        interactionDistance() functions. */
    void setInteractionPoint(int m, double s);

    /** This function stores the specified spatial cell index, distance and cumulative optical
        depth to the initial position of the interaction point for later retrieval through the
        interactionCellIndex(), interactionDistance(), and interactionOpticalDepth() functions. */
    void setInteractionPoint(int m, double s, double tau);

    /** This function returns the spatial cell index corresponding to the interaction point most
        recently calculated by the findInteractionPoint() function or set by the
        setInteractionPoint() function, or -1 if these functions have never been called or if there
        was no interaction point within the path. */
    int interactionCellIndex() const { return _interactionCellIndex; }

    /** This function returns the distance along the path from its initial position to the
        interaction point most recently calculated by the findInteractionPoint() function or set by
        the setInteractionPoint() function, or zero if these functions have never been called or if
        there was no interaction point within the path. */
    double interactionDistance() const { return _interactionDistance; }

    /** This function returns the cumulative optical depth along the path from its initial position
        to the interaction point most recently set by the setInteractionPoint() function, or zero
        if this functions has never been called. */
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
