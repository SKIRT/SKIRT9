/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef TIMEGRID_HPP
#define TIMEGRID_HPP

#include "Range.hpp"
#include "SimulationItem.hpp"

//////////////////////////////////////////////////////////////////////

/** TimeGrid is an abstract class that defines the interface for time grids that can be used, for
    example, to specify the time bins in a light curve instrument.

    A time grid consists of one or more non-overlapping, nonempty time bins in increasing time
    order. Each bin is defined by its left and right borders and has a characteristic time that
    falls inside the bin. The left border is considered to be inside the bin; the right border is
    considered to be outside of the bin. Neighboring bins may have a common border but can also be
    disconnected. Note that time values can be negative.

    Formally, assuming \f$N>0\f$ bins with zero-based indices, we have \f[ t^\mathrm{left}_k \le
    t^\mathrm{c}_k < t^\mathrm{right}_k, \quad k=0\dots N-1 \f] and if \f$N>1\f$, we additionally
    have \f[ t^\mathrm{right}_k \le t^\mathrm{left}_{k+1}, \quad k=0\dots N-2. \f] Finally, each
    bin of course has an associated bin width, \f[t^\mathrm{right}_k - t^\mathrm{left}_k > 0, \quad
    k=0\dots N-1.\f]

    The public interface offers functions to obtain the properties of each bin, and to determine
    the (index of the) bin that contains a given time point.

    A TimeGrid subclass is expected to implement the getTimeBins() function, which is invoked by
    this base class during setup to initialize the time grid. */
class TimeGrid : public SimulationItem
{
    ITEM_ABSTRACT(TimeGrid, SimulationItem, "a time grid")
    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function invokes the getTimeBins() function, implemented in each subclass, to obtain a
        list of time bins. After verifying that the bins are properly formed and sorted, the
        function builds a data structure that allows the bin() function, implemented here, to
        efficiently locate the bin containing a given time point. */
    void setupSelfBefore() override;

    /** This local class represents a single time bin. */
    class Bin
    {
    private:
        double _left, _time, _right;

    public:
        Bin(double left, double time, double right) : _left{left}, _time{time}, _right{right} {}
        double left() const { return _left; }
        double time() const { return _time; }
        double right() const { return _right; }
    };

    /** This function must be implemented in a subclass; it is invoked by this base class during
        setup. The function must place the time bins for this grid in the specified vector, which
        is guaranteed to be empty upon invocation. The bins must be sorted in increasing order and
        must otherwise conform to the restrictions listed in the class header. */
    virtual void getTimeBins(vector<Bin>& bins) const = 0;

    //======================== Public interface =======================

public:
    /** This function returns the number of bins \f$N\f$ in the grid. */
    int numBins() const;

    /** This function returns the characteristic time \f$t^\mathrm{c}_k\f$ for the bin
        corresponding to the index \f$k\f$. */
    double time(int k) const;

    /** This function returns the left border \f$t^\mathrm{left}_k\f$ of the bin corresponding to
        the index \f$k\f$. */
    double left(int k) const;

    /** This function returns the right border \f$t^\mathrm{right}_k\f$ of the bin corresponding to
        the index \f$k\f$. */
    double right(int k) const;

    /** This function returns width \f$t^\mathrm{right}_k - t^\mathrm{left}_k\f$ of the bin
        corresponding to the index \f$k\f$. */
    double width(int k) const;

    /** This function returns the range covered by the time grid, running from the left border of
        the leftmost bin to the right border of the rightmost bin. */
    Range range() const;

    /** This function returns the index \f$k\f$ of the time bin that contains the specified time
        \f$t\f$, i.e. for which \f$t^\mathrm{left}_k \le t < t^\mathrm{right}_k\f$. If no times
        bins match this condition, the function returns -1. */
    int binForTime(double time) const;

    /** This function returns the index \f$k\f$ of the time bin that contains the time
        corresponding to the specified distance \f$d\f$ at the speed of light, i.e. the bin for
        which \f$t^\mathrm{left}_k \le d/c < t^\mathrm{right}_k\f$. If no times bins match this
        condition, the function returns -1. */
    int binForDistance(double distance) const;

    //======================== Data Members ========================

private:
    // initialized in setupSelfBefore()
    vector<Bin> _bins;        // ordered time bins
    vector<double> _borders;  // ordered unique border points (nr of entries depends on whether bins are adjacent)
    vector<int> _indices;     // indices of the  bins defined by the border points, or -1 if not inside bin
};

//////////////////////////////////////////////////////////////////////

#endif
