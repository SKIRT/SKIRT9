/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "LinTimeGrid.hpp"
#include "FatalError.hpp"
#include "NR.hpp"

////////////////////////////////////////////////////////////////////

void LinTimeGrid::getTimeBins(vector<Bin>& bins) const
{
    // verify property values
    if (_maxTime <= _minTime) throw FATALERROR("the maximum time should be larger than the minimum time");

    // construct a grid of characteristic time points
    Array times;
    double halfwidth = 0.5 * NR::buildLinearGrid(times, _minTime, _maxTime, _numTimes - 1);

    // build the list of bins, making sure to use the exact same value for common borders
    bins.reserve(_numTimes);
    bins.emplace_back(times[0] - halfwidth, times[0], times[0] + halfwidth);
    for (int k = 1; k != _numTimes; ++k) bins.emplace_back(bins.back().right(), times[k], times[k] + halfwidth);
}

////////////////////////////////////////////////////////////////////
