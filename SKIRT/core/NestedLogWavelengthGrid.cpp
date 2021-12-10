/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "NestedLogWavelengthGrid.hpp"
#include "FatalError.hpp"
#include "NR.hpp"

////////////////////////////////////////////////////////////////////

void NestedLogWavelengthGrid::setupSelfBefore()
{
    DisjointWavelengthGrid::setupSelfBefore();

    // verify property values
    if (_minWavelengthSubGrid < _minWavelengthBaseGrid || _maxWavelengthSubGrid <= _minWavelengthSubGrid
        || _maxWavelengthBaseGrid < _maxWavelengthSubGrid)
        throw FATALERROR("the high-resolution subgrid should be properly nested in the low-resolution grid");

    // build the high- and low-resolution grids, independently
    Array lambdabasev, lambdasubv;
    NR::buildLogGrid(lambdabasev, _minWavelengthBaseGrid, _maxWavelengthBaseGrid, _numWavelengthsBaseGrid - 1);
    NR::buildLogGrid(lambdasubv, _minWavelengthSubGrid, _maxWavelengthSubGrid, _numWavelengthsSubGrid - 1);

    // merge the two grids (don't worry about order because wavelengths will be sorted later anyway)
    vector<double> lambdav;
    for (double lambda : lambdabasev)
    {
        // use the actual grid limits because they might slightly differ from the specified limits
        if (lambda < lambdasubv[0] || lambda > lambdasubv[_numWavelengthsSubGrid - 1]) lambdav.push_back(lambda);
    }
    for (double lambda : lambdasubv) lambdav.push_back(lambda);

    // store the result
    setWavelengthRange(NR::array(lambdav), true);
}

////////////////////////////////////////////////////////////////////
