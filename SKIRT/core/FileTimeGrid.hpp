/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef FILETIMEGRID_HPP
#define FILETIMEGRID_HPP

#include "DisjointTimeGrid.hpp"

////////////////////////////////////////////////////////////////////

/** FileTimeGrid is a subclass of the DisjointTimeGrid class representing time
    grids loaded from an input file. The floating point numbers in the first column of the text
    file specify the characteristic times. The unit is second.
    In all cases, the times loaded from the file represent the characteristic times of
    the grid bins, and the bin borders are derived automatically. Note that the outermost bin
    borders are placed beyond the outermost characteristic times.

    If \em relativeHalfWidth is zero (the default value) the class constructs a consecutive range
    of adjacent time bins in linear or logarithmic space depending on the value of the \em
    log flag. For more details, refer to the DisjointTimeGrid::settimeRange() function.

    If \em relativeHalfWidth is nonzero, the class constructs a set of distinct nonadjacent
    time bins with the specified relative half bin width. In this case, the value of the \em
    log flag is ignored. For more details, refer to the DisjointTimeGrid::settimeBins()
    function.

    Note that the grid value of the characteristic times allows the negative time lags,
    which can be useful for setting multiple sources. If you would like to use logarithmic spacing for negative time lags,
    you can criate a logarithmic grid which is positive and then shift the grid to negative values by setting the time lags
    to negative values in the input file.
    */
class FileTimeGrid : public DisjointTimeGrid
{
    ITEM_CONCRETE(FileTimeGrid, DisjointTimeGrid, "a time grid loaded from a text file")
        ATTRIBUTE_TYPE_DISPLAYED_IF(FileTimeGrid, "Level2")

        PROPERTY_STRING(filename, "the name of the file with the characteristic time")

        PROPERTY_DOUBLE(relativeHalfWidth, "the relative half width for discrete bins, or zero for a consecutive range")
        ATTRIBUTE_MIN_VALUE(relativeHalfWidth, "[0")
        ATTRIBUTE_MAX_VALUE(relativeHalfWidth, "1[")
        ATTRIBUTE_DEFAULT_VALUE(relativeHalfWidth, "0")

        PROPERTY_BOOL(log, "use logarithmic scale")
        ATTRIBUTE_DEFAULT_VALUE(log, "true")
        ATTRIBUTE_RELEVANT_IF(log, "!relativeHalfWidth")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function reads the time grid points from the specified file. */
    void setupSelfBefore() override;
};

////////////////////////////////////////////////////////////////////

#endif
