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
    WavelengthGrid::setupSelfBefore();

    // verify property values
    if (_minWavelengthSubGrid <= _minWavelengthBaseGrid ||
        _maxWavelengthSubGrid <= _minWavelengthSubGrid ||
        _maxWavelengthBaseGrid <= _maxWavelengthSubGrid)
        throw FATALERROR("the high-resolution subgrid should be properly nested in the low-resolution grid");

    // build the high- and low-resolution grids, independently
    Array lambdalowv, lambdazoomv;
    NR::buildLogGrid(lambdalowv, _minWavelengthBaseGrid, _maxWavelengthBaseGrid, _numWavelengthsBaseGrid-1);
    NR::buildLogGrid(lambdazoomv, _minWavelengthSubGrid, _maxWavelengthSubGrid, _numWavelengthsSubGrid-1);

    // merge the two grids
    vector<double> lambdav;
    for (int ell=0; ell<_numWavelengthsBaseGrid; ell++)
    {
        double lambda = lambdalowv[ell];
        if (lambda<_minWavelengthSubGrid) lambdav.push_back(lambda);
    }
    for (int ell=0; ell<_numWavelengthsSubGrid; ell++)
    {
        double lambda = lambdazoomv[ell];
        lambdav.push_back(lambda);
    }
    for (int ell=0; ell<_numWavelengthsBaseGrid; ell++)
    {
        double lambda = lambdalowv[ell];
        if (lambda>_maxWavelengthSubGrid) lambdav.push_back(lambda);
    }

    // store the result
    setWavelengths(NR::array(lambdav));
}

////////////////////////////////////////////////////////////////////
