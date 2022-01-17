/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "LogWavelengthGrid.hpp"
#include "FatalError.hpp"
#include "NR.hpp"

////////////////////////////////////////////////////////////////////

void LogWavelengthGrid::setupSelfBefore()
{
    DisjointWavelengthGrid::setupSelfBefore();

    // verify property values
    if (_maxWavelength <= _minWavelength) throw FATALERROR("the longest wavelength should be larger than the shortest");

    // construct the grid
    Array lambdav;
    NR::buildLogGrid(lambdav, _minWavelength, _maxWavelength, _numWavelengths - 1);
    setWavelengthRange(lambdav, true);
}

////////////////////////////////////////////////////////////////////
