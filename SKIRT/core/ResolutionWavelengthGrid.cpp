/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "ResolutionWavelengthGrid.hpp"
#include "FatalError.hpp"
#include "NR.hpp"

////////////////////////////////////////////////////////////////////

void ResolutionWavelengthGrid::setupSelfBefore()
{
    DisjointWavelengthGrid::setupSelfBefore();

    // verify property values
    if (_maxWavelength <= _minWavelength) throw FATALERROR("the longest wavelength should be larger than the shortest");

    // determine the number of bins - 1
    int numWavelengths = std::ceil(std::log(_maxWavelength / _minWavelength) / std::log(1. + 1. / _resolution));

    // construct the grid
    Array lambdav;
    NR::buildLogGrid(lambdav, _minWavelength, _maxWavelength, numWavelengths);
    setWavelengthRange(lambdav, true);
}

////////////////////////////////////////////////////////////////////
