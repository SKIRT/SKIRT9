/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "LinWavelengthGrid.hpp"
#include "FatalError.hpp"
#include "NR.hpp"

////////////////////////////////////////////////////////////////////

void LinWavelengthGrid::setupSelfBefore()
{
    DisjointWavelengthGrid::setupSelfBefore();

    // verify property values
    if (_maxWavelength <= _minWavelength) throw FATALERROR("the longest wavelength should be larger than the shortest");

    // construct the grid
    Array lambdav;
    NR::buildLinearGrid(lambdav, _minWavelength, _maxWavelength, _numWavelengths - 1);
    setWavelengthRange(lambdav, false);
}

////////////////////////////////////////////////////////////////////
