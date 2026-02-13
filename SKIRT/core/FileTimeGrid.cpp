/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "FileTimeGrid.hpp"
#include "TextInFile.hpp"

////////////////////////////////////////////////////////////////////

void FileTimeGrid::setupSelfBefore()
{
    DisjointTimeGrid::setupSelfBefore();

    // read the times from the input file
    TextInFile infile(this, _filename, "time grid");
    infile.addColumn("time", "time", "s");
    Array timeLags;
    infile.readAllColumns(timeLags);
    infile.close();

    setTimeRange(timeLags, _log);
}

//////////////////////////////////////////////////////////////////////
