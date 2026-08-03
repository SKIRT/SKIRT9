/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "LogTimeGrid.hpp"
#include "FatalError.hpp"
#include "NR.hpp"

////////////////////////////////////////////////////////////////////

void LogTimeGrid::getTimeBins(vector<Bin>& bins) const
{
    // verify property values
    if (_maxTime <= _minTime) throw FATALERROR("the maximum time should be larger than the minimum time");

    // construct a grid of characteristic time points
    Array times;
    NR::buildLogGrid(times, _minTime, _maxTime, _numTimes - 1);

    // build the list of bins, making sure to use the exact same value for common borders
    bins.reserve(_numTimes);
    {
        double left = times[0] * sqrt(times[0] / times[1]);
        double right = sqrt(times[1] * times[0]);
        bins.emplace_back(left - _offset, times[0] - _offset, right - _offset);
    }
    for (int k = 1; k != _numTimes - 1; ++k)
    {
        double right = sqrt(times[k + 1] * times[k]);
        bins.emplace_back(bins.back().right(), times[k] - _offset, right - _offset);
    }
    {
        double right = times[_numTimes - 1] * sqrt(times[_numTimes - 1] / times[_numTimes - 2]);
        bins.emplace_back(bins.back().right(), times[_numTimes - 1] - _offset, right - _offset);
    }
}

////////////////////////////////////////////////////////////////////
