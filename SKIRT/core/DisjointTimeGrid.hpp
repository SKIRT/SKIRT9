/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef DISJOINTTIMEGRID_HPP
#define DISJOINTTIMEGRID_HPP

#include "TimeGrid.hpp"

//////////////////////////////////////////////////////////////////////

/** DisjointTimeGrid is an abstract class that represents time grids with
    straightforward, non-overlapping bins with constant transmission across each bin.

    Specifically, a disjoint time grid consists of non-overlapping but possibly adjacent
    time bins in increasing time order, with constant maximum transmission within the
    bins and zero transmission outside of the bins. Each bin is defined by its left and right
    borders and has a characteristic time that falls inside the bin. The left border is
    considered to be inside of the bin; the right border is considered to be outside of the bin.
    Neighboring bins may have a common border but can also be disconnected.

    Formally, assuming \f$N>0\f$ bins with zero-based indices, we have \f[
    \time^\mathrm{left}_\ell \le \time^\mathrm{c}_\ell < \time^\mathrm{right}_\ell, \quad
    \ell=0\dots N-1 \f] and if \f$N>1\f$, we additionally have \f[ \time^\mathrm{right}_\ell \le
    \time^\mathrm{left}_{\ell+1}, \quad \ell=0\dots N-2. \f] Finally, each bin of course has an
    associated bin width, \f[\time^\mathrm{right}_\ell - \time^\mathrm{left}_\ell > 0, \quad
    \ell=0\dots N-1.\f]

    A DisjointTimeGrid subclass is expected to invoke one of the settimeXXX() functions
    during setup to initialize the time grid. The current implementation offers two such
    functions: one to specify a consecutive range of adjacent time bins given a list of
    characteric times, and another one to specify distinct, nonadjacent time bins given
    a list of characteric times and a relative bin width. Other options can be added as the
    need arises. */
class DisjointTimeGrid : public TimeGrid
{
    ITEM_ABSTRACT(DisjointTimeGrid, TimeGrid, "a time grid with non-overlapping bins")
    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function verifies that the time bins have been initialized by a subclass calling
        one of the settimeXXX() functions of this class in their setupSelfBefore() function.
        */
    void setupSelfAfter() override;

    /** This function initializes the time grid to a consecutive range of adjacent time
        bins given a list of characteric times. The subclass determines a list of characteric
        times according to some predefined scheme, and the bin borders and bin widths are
        automatically determined from that list by this function. If the specified time list
        is empty, or if there are duplicate values (which would lead to empty bins), the function
        throws a fatal error.

        The function first sorts the specified characteristic times in ascending order and
        then calculates the bin borders assuming linear or logarithmic scaling depending on the
        value of the \em logScale flag. The inner border between two consecutive characteristic
        times is placed at the mid-point of those two times (in linear or logarithmic
        space), and the outer borders at the edges of the time range are placed such that the
        outer bins have the same width as the respective adjacent bins (in linear or logarithmic
        space). If the time grid has just a single time, the borders are placed just
        next to the time to form a narrow bin. Finally, the function trivially calculates the
        time bin widths from the bin borders.

        For linear scaling, the corresponding formulas are trivial. For logarithmic scaling, the
        formalas in logarithmic space translate easily to equivalent but more efficient formulas in
        real space. For the inner borders this yields the geometric mean of the two adjacent
        characteristic times, i.e. \f$\time^\mathrm{right}_{\ell-1} =
        \time^\mathrm{left}_\ell = \sqrt{\time^\mathrm{c}_{\ell-1}\time^\mathrm{c}_\ell}\;,
        \ell=1\dots N-1\f$. The leftmost outer border is placed at \f$\time^\mathrm{left}_0 =
        \sqrt{(\time^\mathrm{c}_{0})^3/\time^\mathrm{c}_1}\f$, and the rightmost outer border
        is placed at \f$\time^\mathrm{right}_{N-1} =
        \sqrt{(\time^\mathrm{c}_{N-1})^3/\time^\mathrm{c}_{N-2}}\f$.

        If there is just a single time in the grid, the outer borders are placed (for both
        linear and logarithmic scaling) according to \f$\time^\mathrm{left}_0 =
        \time^\mathrm{c}_{0}(1-1/1000)\f$ and \f$\time^\mathrm{right}_0 =
        \time^\mathrm{c}_{0}(1+1/1000)\f$. */
    void setTimeRange(const Array& timev, bool logScale);

    /** This function initializes the time grid to a set of distinct, nonadjacent time
        bins given a list of characteric times and a relative half bin width. The subclass
        determines a list of characteric times and a relative half bin width, and the bin
        borders and bin widths are automatically calculated from that information by this function.
        If the specified time list is empty, or if the calculated bins overlap, the function throws a fatal error.

        Specifically, the function first sorts the specified characteristic times in
        ascending order and then calculates the bin borders using \f$\time^\mathrm{left}_\ell =
        \time^\mathrm{c}_\ell(1-w)\f$ and \f$\time^\mathrm{right}_\ell =
        \time^\mathrm{c}_\ell(1+w)\;, \ell=0\dots N-1\f$, where \f$w\f$ is the specified relative
        half bin width. If the \em constantWidth flag is true, the width for the shortest
        time is used for all bin widths instead. Finally the function trivially calculates
        the time bin widths from the bin borders. */
    void setTimeBins(const Array& timev, double relativeHalfWidth, bool constantWidth = false);

    /** This function initializes the time grid to a consecutive range of \f$N>0\f$ adjacent
        time bins given a list of \f$N+1\f$ time bin borders. The subclass determines a
        list of bin borders according to some predefined scheme, and the characteristic times
        and bin widths are automatically determined from that list by this function. If the
        specified list has fewer than two bin borders, or if there are duplicate values (which
        would lead to empty bins), the function throws a fatal error.

        The function first sorts the specified time bin borders in ascending order and then
        calculates the characteristic times assuming linear scaling (arithmetic mean) or
        logarithmic scaling (geometric mean) depending on the value of the \em logScale flag. */
    void setTimeBorders(const Array& borderv, bool logScale);

    /** This function initializes the time grid from a list of interleaved bin border points
        and corresponding characteristic times. (i.e., borders and characteristic times
        alternate). The number of values must be uneven and at least three. The list must be in
        strictly increasing or decreasing order, which means duplicates are not allowed, except
        that a zero characteristic time indicates a segment that is not part of the grid,
        i.e. that lies between two non-adjacent bins. In other words, this option allows to (1)
        arbitrarily place characteristic times within each bin and (2) to specify
        intermediate time ranges that are not covered by any bin. */
    void setTimeSegments(const Array& bordcharv);

    /** This function initializes the time grid from a list of bin border points and a
        corresponding list of characteristic times. A zero characteristic time
        indicates a segment that is not part of the grid, i.e. that lies between two non-adjacent
        bins or beyond the outer grid borders. In other words, the subclass has full control over
        the placement of bin borders and characteristic times.

        The two lists must have the same size. The borders must be listed in strictly increasing
        order of time (i.e. there cannot be any duplicates), all nonzero characteristic
        times must lie within their corresponding bin, the last characteristic time
        must be zero, and there must be at least one nonzero characteristic time. If these
        requirements are violated, the behavior of this function is undefined. */
    void setTimeSegments(const vector<double>& borderv, const vector<double>& characv);

    //================= Functions implementing virtual base class functions ===================

public:
    /** This function returns the number of bins, \f$N\f$, in the grid (or equivalently, the number
        of characteristic times). */
    int numBins() const override;

    /** This function returns the characteristic time \f$\time^\mathrm{c}_\ell\f$
       corresponding to the index \f$\ell\f$. */
    double time(int ell) const override;

    /** This function returns the left border of the time bin corresponding to the index
        \f$\ell\f$, i.e. \f$\time^\mathrm{left}_\ell\f$. */
    double leftBorder(int ell) const override;

    /** This function returns the right border of the time bin corresponding to the index
        \f$\ell\f$, i.e. \f$\time^\mathrm{right}_\ell\f$. */
    double rightBorder(int ell) const override;

    /** This function returns the width of the time bin corresponding to the index
        \f$\ell\f$, i.e. \f$\time^\mathrm{right}_\ell - \time^\mathrm{left}_\ell\f$. */
    double effectiveWidth(int ell) const override;

    /** This function returns the relative transmission for the time bin corresponding to the
        index \f$\ell\f$ at the time \f$\time\f$. For the present class, it always returns
        1. */
    double transmission(int ell, double time) const override;

    /** This function returns a single-element list with the index \f$\ell\f$ of the time bin
        that contains the specified time \f$\time\f$, i.e. for which
        \f$\time^\mathrm{left}_\ell <= \time < \time^\mathrm{right}_\ell\f$. If \f$\time\f$
        does not lie inside one of the time bins, the function returns an empty list. */
    vector<int> bins(double time) const override;

    /** This function returns the index \f$\ell\f$ of the time bin that contains the
        specified time \f$\time\f$, i.e. for which \f$\time^\mathrm{left}_\ell <= \time
        < \time^\mathrm{right}_\ell\f$. If \f$\time\f$ does not lie inside one of the
        time bins, the function returns -1. */
    int bin(double time) const override;

    //=============== Functions specific to disjoint time grids =================

public:
    /** This function returns (a reference to) the list of characteristic times in this
        time grid. In combination with the dtimev() function, it allows easily expressing
        calculations involving consecutive time grids. */
    const Array& timev() const { return _timev; }

    /** This function returns (a reference to) the list of bin widths in this time grid. In
        combination with the timev() function, it allows easily expressing calculations involving
        consecutive time grids. */
    const Array& dtimev() const { return _dtimev; }

    /** This function returns a list of the characteristic times in this time grid
        extended with the outermost bin border point on each side. The list has thus two additional
        points, one on each side, and as a result covers the complete time range of the grid,
        including the widths of the outer bins. This extended time list can be used in
        situations where one needs to calculate/interpolate some function over the complete range
        of the time grid and not just up to the outermost characteristic times. */
    Array exttimev() const;

    /** This function returns a list of bin widths in this time grid extended with a zero
        value on each side. This extended time bin width list can be used, for example, to
        integrate a function discretized on the extended time grid returned by the
        exttimev() function over the time range. */
    Array extdtimev() const;

    //======================== Data Members ========================

private:
    // subclasses should call settimeXXX() in their setupSelfBefore() to initialize these tables
    Array _timev;       // N characteristic times
    Array _dtimev;      // N time bin widths
    Array _timeleftv;   // N left time bin borders
    Array _timerightv;  // N right time bin widths
    Array _borderv;       // K=N+1 or K=N*2 ordered border points (depending on whether bins are adjacent)
    vector<int> _ellv;    // K+1 indices of the time bins defined by the border points, or -1 if out of range
};

//////////////////////////////////////////////////////////////////////

#endif
