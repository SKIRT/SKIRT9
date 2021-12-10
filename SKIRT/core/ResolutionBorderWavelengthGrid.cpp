/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "ResolutionBorderWavelengthGrid.hpp"
#include "FatalError.hpp"
#include "NR.hpp"

////////////////////////////////////////////////////////////////////

void ResolutionBorderWavelengthGrid::setupSelfBefore()
{
    DisjointWavelengthGrid::setupSelfBefore();

    // verify property values
    if (_maxWavelength <= _minWavelength) throw FATALERROR("the longest wavelength should be larger than the shortest");

    // determine the number of bins
    int numWavelengthBins = std::ceil(std::log(_maxWavelength / _minWavelength) / std::log(1. + 1. / _resolution));

    // construct the grid
    Array borderv;
    NR::buildLogGrid(borderv, _minWavelength, _maxWavelength, numWavelengthBins);
    setWavelengthBorders(borderv, true);
}

////////////////////////////////////////////////////////////////////
