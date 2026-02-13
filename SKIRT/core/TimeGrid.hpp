/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef TIMEGRID_HPP
#define TIMEGRID_HPP

#include "Array.hpp"
#include "Range.hpp"
#include "SimulationItem.hpp"

//////////////////////////////////////////////////////////////////////

/** TimeGrid is an abstract class that defines the interface for time grids that can be
    used, for example, to specify the time bins used for detecting photon packets in
    instruments.

    A time grid consists of \f$N>0\f$ possibly overlapping time bins. Generally
    speaking, each of these bins is defined through a transmission determining the
    contribution fraction of a detected photon packet to the bin as a function of time. This
    class defines the public interface for a time grid in these general terms. Subclasses can
    implement several bin types, including bins that mimic a particular broadband, or
    straightforward bins with constant transmission across some time interval.

    Key properties of a time bin include its left and right borders, defining a time
    interval, and time lag. The TimeGrid class requires that the
    time lag of a bin falls inside the bin, i.e. \f$t^\mathrm{left}_\ell \le
    t^\mathrm{c}_\ell \le t^\mathrm{right}_\ell, \ell=0\dots N-1\f$, and that no two bins
    have the same time lag. This allows the TimeGrid class to meaningfully
    sort bins in increasing order of time lag, and assign bin indices accordingly.

    The public interface also includes functions to obtain the relative transmission for a given
    bin as a function of time, defined as the transmission at that time divided by the
    maximum transmission value, and to obtain a bin's effective width, defined as the horizontal size of
    a rectangle with height equal to the maximum transmission value and with the same area as the one
    covered by the band's transmission.

    Finally, and most importantly, the public interface offers a function to determine the (indices
    of) the bin(s) that may have a nonzero transmission at a given time. */
class TimeGrid : public SimulationItem
{
    ITEM_ABSTRACT(TimeGrid, SimulationItem, "a time grid")
    ITEM_END()

    //======================== Public interface =======================

public:
    /** This function returns the number of bins \f$N\f$ in the grid, or equivalently, the number
        of time lags. Bins are always sorted in order of increasing time lag. */
    virtual int numBins() const = 0;

    /** This function returns the time lag \f$t^\mathrm{c}_\ell\f$
       corresponding to the index \f$\ell\f$. The time lag of a bin is always
       inside the bin, i.e. \f$t^\mathrm{left}_\ell \le t^\mathrm{c}_\ell \le
       t^\mathrm{right}_\ell\f$. Bins are always sorted in order of
       increasing time lag. */
    virtual double time(int ell) const = 0;

    /** This function returns the left border of the time bin corresponding to the index
        \f$\ell\f$, i.e. \f$t^\mathrm{left}_\ell\f$. The transmission for this bin is
        guaranteed to be zero for all times shorter than the left border. */
    virtual double leftBorder(int ell) const = 0;

    /** This function returns the right border of the time bin corresponding to the index
        \f$\ell\f$, i.e. \f$t^\mathrm{right}_\ell\f$. The transmission for this bin is
        guaranteed to be zero for all times longer than the right border. */
    virtual double rightBorder(int ell) const = 0;

    /** This function returns the effective width of the time bin corresponding to the index
        \f$\ell\f$. The effective width is defined as the horizontal size of a rectangle with
        height equal to the maximum transmission and with the same area as the one covered by the
        band's transmission curve. For bins with a constant transmission over the complete
        interval, the effective width is simply \f$t^\mathrm{right}_\ell -
        t^\mathrm{left}_\ell\f$. For bins with non-constant transmission, this is no longer
        true. */
    virtual double effectiveWidth(int ell) const = 0;

    /** This function returns the relative transmission for the time bin corresponding to the
        index \f$\ell\f$ at the time \f$t\f$. The relative transmission is defined as
        the transmission at that time divided by the maximum transmission for the bin. */
    virtual double transmission(int ell, double time) const = 0;

    /** This function returns a list of indices \f$\ell_k\f$ of the time bins that may have a
        nonzero transmission at the specified time \f$t\f$, i.e. for which
        \f$t^\mathrm{left}_\ell \le t \le t^\mathrm{right}_\ell\f$. If no
        times bins match this condition, the function returns an empty list. */
    virtual vector<int> bins(double time) const = 0;

    /** This function returns the index \f$\ell\f$ of one the time bins that may have a
        nonzero transmission at the specified time \f$t\f$, i.e. for which
        \f$t^\mathrm{left}_\ell \le t \le t^\mathrm{right}_\ell\f$. If no
        times bins match this condition, the function returns -1. If multiple bins match this
        condition, the function returns the index for the bin with the shortest characteristic
        time. */
    virtual int bin(double time) const = 0;

    /** This function returns the time range covered by the time grid, which is defined
        as the range from the left border of the leftmost bin to the right border of the rightmost
        bin. This range includes all times possibly covered by the time grid except in
        the rare case of overlapping bins where an inner bin is so wide that it outer limit extends
        beyond the outer bin. */
    Range timeRange() const;
};

//////////////////////////////////////////////////////////////////////

#endif
