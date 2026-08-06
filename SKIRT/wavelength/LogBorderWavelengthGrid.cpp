/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "LogBorderWavelengthGrid.hpp"
#include "FatalError.hpp"
#include "NR.hpp"

////////////////////////////////////////////////////////////////////

void LogBorderWavelengthGrid::setupSelfBefore()
{
    DisjointWavelengthGrid::setupSelfBefore();

    // verify property values
    if (_maxWavelength <= _minWavelength) throw FATALERROR("the longest wavelength should be larger than the shortest");

    // construct the grid
    Array borderv;
    NR::buildLogGrid(borderv, _minWavelength, _maxWavelength, _numWavelengthBins);
    setWavelengthBorders(borderv, true);
}

////////////////////////////////////////////////////////////////////
