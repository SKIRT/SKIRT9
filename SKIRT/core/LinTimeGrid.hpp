/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef LINTIMEGRID_HPP
#define LINTIMEGRID_HPP

#include "TimeGrid.hpp"

//////////////////////////////////////////////////////////////////////

/** LinTimeGrid represents linearly distributed time grids. The characteristic times of the grid
    bins are equally distributed between and including the specified minimum and maximum time.
    These minimum and/or maximum time values are allowed to be negative.

    The outermost bins are given the same width as the inner bins, which implies that the outermost
    bin borders are placed half a bin width beyond the specified minimum and maximum time. The grid
    must have at least two bins, which then have the specified minimum and maximum time as their
    respective characteristic time. */
class LinTimeGrid : public TimeGrid
{
    ITEM_CONCRETE(LinTimeGrid, TimeGrid, "a linear time grid")
        ATTRIBUTE_TYPE_DISPLAYED_IF(FileTimeGrid, "Level2")

        PROPERTY_DOUBLE(minTime, "the minimum time")
        ATTRIBUTE_QUANTITY(minTime, "timelag")

        PROPERTY_DOUBLE(maxTime, "the maximum time")
        ATTRIBUTE_QUANTITY(maxTime, "timelag")

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
