/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef FILETIMEGRID_HPP
#define FILETIMEGRID_HPP

#include "TimeGrid.hpp"

////////////////////////////////////////////////////////////////////

/** FileTimeGrid represents time grids loaded from an input file. The floating point numbers in the
    first three columns of the text file specify respectively the characteristic time, the left
    border, and the right border of each bin. The default unit is second (s), but this can be
    overridden by a column header (see TextInFile). Any additional columns in the file are ignored.

    The bins must be non-empty and non-overlapping, and must be sorted in increasing time order. It
    is allowed to have gaps between the bins. For a formal statement of the requirements, see
    TimeGrid. Note that time values are allowed to be negative.

    \note When specifying consecutive "touching" bins, make sure that the right border of the first
    bin is exactly equal to the left border of the second bin to the precision of the values listed
    in the input file. if not, the import might fail because the bins are considered to overlap, or
    there will be an unintended (small) gap in the time grid.

    */
class FileTimeGrid : public TimeGrid
{
    ITEM_CONCRETE(FileTimeGrid, TimeGrid, "a time grid loaded from a text file")
        ATTRIBUTE_TYPE_DISPLAYED_IF(FileTimeGrid, "Level2")

        PROPERTY_STRING(filename, "the name of the file with the time bins")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function is invoked during setup. It places the time bins for this grid in the
        specified vector, which is guaranteed to be empty upon invocation. */
    void getTimeBins(vector<Bin>& bins) const override;
};

////////////////////////////////////////////////////////////////////

#endif
