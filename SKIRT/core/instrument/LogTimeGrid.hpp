/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef LOGTIMEGRID_HPP
#define LOGTIMEGRID_HPP

#include "TimeGrid.hpp"

//////////////////////////////////////////////////////////////////////

/** LogTimeGrid represents logarithmically distributed time grids. The characteristic times of the
    grid bins are equally distributed (in log space) between and including the specified minimum
    and maximum time. These minimum and/or maximum time values are allowed to be negative.

    The outermost bins are given the same width as the inner bins (in log space), which implies
    that the outermost bin borders are placed half a bin width beyond the specified minimum and
    maximum time. The grid must have at least two bins, which then have the specified minimum and
    maximum time as their respective characteristic time.

    Because a logarithmic grid cannot have negative or zero values, the additional \em offset
    property specifies a value that is subtracted from each of the grid points. A positive value
    shifts the complete grid towards and possibly beyond the origin. A negative value shifts
    the grid away from the origin. */
class LogTimeGrid : public TimeGrid
{
    ITEM_CONCRETE(LogTimeGrid, TimeGrid, "a logarithmic time grid")
        ATTRIBUTE_TYPE_DISPLAYED_IF(FileTimeGrid, "Level2")

        PROPERTY_DOUBLE(minTime, "the minimum time")
        ATTRIBUTE_QUANTITY(minTime, "timelag")
        ATTRIBUTE_MIN_VALUE(minTime, "]0")

        PROPERTY_DOUBLE(maxTime, "the maximum time")
        ATTRIBUTE_QUANTITY(maxTime, "timelag")
        ATTRIBUTE_MIN_VALUE(maxTime, "]0")

        PROPERTY_DOUBLE(offset, "the global time offset subtracted from the grid")
        ATTRIBUTE_QUANTITY(offset, "timelag")
        ATTRIBUTE_DEFAULT_VALUE(offset, "0")

        PROPERTY_INT(numTimes, "the number of time grid points")
        ATTRIBUTE_MIN_VALUE(numTimes, "2")
        ATTRIBUTE_DEFAULT_VALUE(numTimes, "25")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function is invoked during setup. It places the time bins for this grid in the
        specified vector, which is guaranteed to be empty upon invocation. */
    void getTimeBins(vector<Bin>& bins) const override;
};

//////////////////////////////////////////////////////////////////////

#endif
