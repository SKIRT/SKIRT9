/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "TimeGrid.hpp"
#include "Constants.hpp"
#include "FatalError.hpp"
#include "NR.hpp"

//////////////////////////////////////////////////////////////////////

void TimeGrid::setupSelfBefore()
{
    SimulationItem::setupSelfBefore();

    // get time bins from subclass
    getTimeBins(_bins);

    // --- verify conformance ---

    // verify that there is at least one bin
    int numBins = _bins.size();
    if (numBins == 0) throw FATALERROR("There must be at least on time bin");

    // verify that the bin borders and characteristic time are properly ordered
    for (const auto& bin : _bins)
    {
        if (bin.left() > bin.time() || bin.time() >= bin.right())
            throw FATALERROR("Characteristic time must be between left and right bin borders");
    }

    // verify that the bins do not overlap and are in increasing time order
    for (int k = 1; k != numBins; ++k)
    {
        if (_bins[k - 1].right() > _bins[k].left())
            throw FATALERROR("Time bins overlap or are not in increasing time order");
    }

    // --- build structure for binary search ---
    // build a list of unique borders and corresponding bin indices or -1 for gaps between bins

    // always simply add the first bin
    _borders.push_back(_bins[0].left());
    _borders.push_back(_bins[0].right());
    _indices.push_back(0);

    // loop over remaining bins, if any
    for (int k = 1; k != numBins; ++k)
    {
        const auto& bin = _bins[k];
        if (_borders.back() == bin.left())
        {
            // this bin is adjacent to the previous one
            _borders.push_back(bin.right());
            _indices.push_back(k);
        }
        else
        {
            // there is a gap between this bin and the previous one
            _indices.push_back(-1);
            _borders.push_back(bin.left());
            _borders.push_back(bin.right());
            _indices.push_back(k);
        }
    }

    // add another "gap" for the space beyond the last bin
    _indices.push_back(-1);
}

//////////////////////////////////////////////////////////////////////

int TimeGrid::numBins() const
{
    return _bins.size();
}

//////////////////////////////////////////////////////////////////////

double TimeGrid::time(int k) const
{
    return _bins[k].time();
}

//////////////////////////////////////////////////////////////////////

double TimeGrid::left(int k) const
{
    return _bins[k].left();
}

//////////////////////////////////////////////////////////////////////

double TimeGrid::right(int k) const
{
    return _bins[k].right();
}

//////////////////////////////////////////////////////////////////////

double TimeGrid::width(int k) const
{
    return _bins[k].right() - _bins[k].left();
}

//////////////////////////////////////////////////////////////////////

Range TimeGrid::range() const
{
    return Range(_bins[0].left(), _bins[_bins.size() - 1].right());
}

//////////////////////////////////////////////////////////////////////

int TimeGrid::binForTime(double time) const
{
    int i = NR::locate(_borders, time);
    return i >= 0 ? _indices[i] : -1;
}

//////////////////////////////////////////////////////////////////////

int TimeGrid::binForDistance(double distance) const
{
    return binForTime(distance / Constants::c());
}

//////////////////////////////////////////////////////////////////////
