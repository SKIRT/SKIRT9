/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "DisjointTimeGrid.hpp"
#include "FatalError.hpp"
#include "NR.hpp"

////////////////////////////////////////////////////////////////////

void DisjointTimeGrid::setupSelfAfter()
{
    TimeGrid::setupSelfAfter();

    // the subclass should have added time
    if (!_timev.size()) throw FATALERROR("Time grid should have been initialized by subclass");
}

////////////////////////////////////////////////////////////////////

void DisjointTimeGrid::setTimeRange(const Array& timev, bool logScale)
{
    // copy and sort the specified characteristic time
    _timev = timev;
    std::sort(begin(_timev), end(_timev));
    size_t n = _timev.size();

    // verify that there is at least one time
    if (!n) throw FATALERROR("There must be at least one time in the grid");

    // verify that there are no duplicates
    if (std::unique(begin(_timev), end(_timev)) != end(_timev))
        throw FATALERROR("There should be no duplicate time in the grid");

    // calculate the bin borders
    _timeleftv.resize(n);
    _timerightv.resize(n);
    _borderv.resize(n + 1);

    if (n == 1)
    {
        // just a single time -> form narrow bin
        _timeleftv[0] = _borderv[0] = _timev[0] * 0.999;
        _timerightv[0] = _borderv[1] = _timev[0] * 1.001;
    }
    else
    {
        if (logScale)
        {
            _timeleftv[0] = _borderv[0] = sqrt(_timev[0] * _timev[0] * _timev[0] / _timev[1]);
            for (size_t ell = 1; ell != n; ++ell)
            {
                _timerightv[ell - 1] = _timeleftv[ell] = _borderv[ell] = sqrt(_timev[ell - 1] * _timev[ell]);
            }
            _timerightv[n - 1] = _borderv[n] =
                sqrt(_timev[n - 1] * _timev[n - 1] * _timev[n - 1] / _timev[n - 2]);
        }
        else
        {
            _timeleftv[0] = _borderv[0] = (3. * _timev[0] - _timev[1]) / 2.;
            for (size_t ell = 1; ell != n; ++ell)
            {
                _timerightv[ell - 1] = _timeleftv[ell] = _borderv[ell] = (_timev[ell - 1] + _timev[ell]) / 2.;
            }
            _timerightv[n - 1] = _borderv[n] = (3. * _timev[n - 1] - _timev[n - 2]) / 2.;
        }
    }

    // calculate the bin widths
    _dtimev = _timerightv - _timeleftv;

    // setup the mapping from border bin indices to actual time bin indices (see bin() function)
    _ellv.resize(n + 2);
    _ellv[0] = -1;
    for (size_t ell = 0; ell != n; ++ell) _ellv[ell + 1] = ell;
    _ellv[n + 1] = -1;
}

////////////////////////////////////////////////////////////////////

void DisjointTimeGrid::setTimeBins(const Array& timev, double relativeHalfWidth, bool constantWidth)
{
    // copy and sort the specified characteristic time
    _timev = timev;
    std::sort(begin(_timev), end(_timev));
    size_t n = _timev.size();

    // verify that there is at least one time
    if (!n) throw FATALERROR("There must be at least one time in the grid");

    // calculate the bin borders
    _timeleftv.resize(n);
    _timerightv.resize(n);
    _borderv.resize(2 * n);
    if (!constantWidth)
    {
        for (size_t ell = 0; ell != n; ++ell)
        {
            _borderv[2 * ell] = _timeleftv[ell] = _timev[ell] * (1. - relativeHalfWidth);
            _borderv[2 * ell + 1] = _timerightv[ell] = _timev[ell] * (1. + relativeHalfWidth);
        }
    }
    else
    {
        double delta = _timev[0] * relativeHalfWidth;
        for (size_t ell = 0; ell != n; ++ell)
        {
            _borderv[2 * ell] = _timeleftv[ell] = _timev[ell] - delta;
            _borderv[2 * ell + 1] = _timerightv[ell] = _timev[ell] + delta;
        }
    }

    // verify that the bins do not overlap
    if (!std::is_sorted(begin(_borderv), end(_borderv)))
        throw FATALERROR("Non-adjacent time bins should not overlap");

    // calculate the bin widths
    _dtimev = _timerightv - _timeleftv;

    // setup the mapping from border bin indices to actual time bin indices (see bin() function)
    _ellv.resize(2 * n + 1);
    _ellv[0] = -1;
    for (size_t ell = 0; ell != n; ++ell)
    {
        _ellv[2 * ell + 1] = ell;
        _ellv[2 * ell + 2] = -1;  // regions between the bins are considered out of range
    }
}

////////////////////////////////////////////////////////////////////

void DisjointTimeGrid::setTimeBorders(const Array& borderv, bool logScale)
{
    // copy and sort the specified bin borders
    _borderv = borderv;
    std::sort(begin(_borderv), end(_borderv));

    // verify that there are at least two bin borders
    if (_borderv.size() < 2) throw FATALERROR("There must be at least two time bin borders in the grid");

    // verify that there are no duplicates
    if (std::unique(begin(_borderv), end(_borderv)) != end(_borderv))
        throw FATALERROR("There should be no duplicate time bin borders in the grid");

    // make n refer to the number of bins, not the number of borders
    size_t n = _borderv.size() - 1;

    // copy the left and right bin borders
    _timeleftv.resize(n);
    _timerightv.resize(n);
    for (size_t ell = 0; ell != n; ++ell)
    {
        _timeleftv[ell] = _borderv[ell];
        _timerightv[ell] = _borderv[ell + 1];
    }

    // calculate the characteristic time
    _timev.resize(n);
    if (logScale)
    {
        for (size_t ell = 0; ell != n; ++ell)
        {
            _timev[ell] = sqrt(_timeleftv[ell] * _timerightv[ell]);
        }
    }
    else
    {
        for (size_t ell = 0; ell != n; ++ell)
        {
            _timev[ell] = (_timeleftv[ell] + _timerightv[ell]) / 2.;
        }
    }

    // calculate the bin widths
    _dtimev = _timerightv - _timeleftv;

    // setup the mapping from border bin indices to actual time bin indices (see bin() function)
    _ellv.resize(n + 2);
    _ellv[0] = -1;
    for (size_t ell = 0; ell != n; ++ell) _ellv[ell + 1] = ell;
    _ellv[n + 1] = -1;
}

////////////////////////////////////////////////////////////////////

void DisjointTimeGrid::setTimeSegments(const Array& bordcharv)
{
    // verify that the number of values is uneven and least three
    size_t n = bordcharv.size();
    if (n < 3) throw FATALERROR("There must be at least three time values in the list");
    if (n % 2 == 0) throw FATALERROR("There number of time values in the list must be uneven");

    // copy the values into vectors that can be passed to the function that will do the actual work
    vector<double> borderv;
    vector<double> characv;
    size_t i = 0;
    while (true)
    {
        borderv.push_back(bordcharv[i++]);
        if (i == n) break;
        characv.push_back(bordcharv[i++]);
    }

    // reverse the lists if required
    if (bordcharv[0] > bordcharv[n - 1])
    {
        std::reverse(borderv.begin(), borderv.end());
        std::reverse(characv.begin(), characv.end());
    }

    // add the "characteristic time" for the final segment outside the grid
    characv.push_back(0.);

    // verify the ordering
    n = borderv.size() - 1;
    for (size_t i = 0; i != n; ++i)
    {
        if (borderv[i + 1] <= borderv[i])
            throw FATALERROR("Time bin borders must be in strictly monotonous order");
        if (std::isinf(characv[i])) characv[i] = 0.;  // handle zero frequencies
        if (characv[i] != 0. && (characv[i] <= borderv[i] || characv[i] >= borderv[i + 1]))
            throw FATALERROR("Characteristic time must be within bin borders");
    }

    // call the function that will do the actual work
    setTimeSegments(borderv, characv);
}

////////////////////////////////////////////////////////////////////

void DisjointTimeGrid::setTimeSegments(const vector<double>& borderv, const vector<double>& characv)
{
    // determine the number of borders and the number of actual bins
    int numBorders = borderv.size();
    int numBins = std::count_if(characv.begin(), characv.end(), [](double c) { return c > 0.; });

    // allocate memory for all arrays
    _timev.resize(numBins);       // index ell
    _dtimev.resize(numBins);      // index ell
    _timeleftv.resize(numBins);   // index ell
    _timerightv.resize(numBins);  // index ell
    _borderv.resize(numBorders);    // index k
    _ellv.resize(numBorders + 1);   // index k+1

    // copy the data into our arrays, and construct the bin index mapping
    _ellv[0] = -1;
    int ell = 0;  // ell is bin index; k is border index
    for (int k = 0; k != numBorders; ++k)
    {
        _borderv[k] = borderv[k];
        if (characv[k] > 0.)
        {
            _timev[ell] = characv[k];
            _dtimev[ell] = borderv[k + 1] - borderv[k];
            _timeleftv[ell] = borderv[k];
            _timerightv[ell] = borderv[k + 1];
            _ellv[k + 1] = ell++;
        }
        else
        {
            _ellv[k + 1] = -1;
        }
    }
}

////////////////////////////////////////////////////////////////////

int DisjointTimeGrid::numBins() const
{
    return _timev.size();
}

////////////////////////////////////////////////////////////////////

double DisjointTimeGrid::time(int ell) const
{
    return _timev[ell];
}

////////////////////////////////////////////////////////////////////

double DisjointTimeGrid::leftBorder(int ell) const
{
    return _timeleftv[ell];
}

////////////////////////////////////////////////////////////////////

double DisjointTimeGrid::rightBorder(int ell) const
{
    return _timerightv[ell];
}

////////////////////////////////////////////////////////////////////

double DisjointTimeGrid::effectiveWidth(int ell) const
{
    return _dtimev[ell];
}

////////////////////////////////////////////////////////////////////

double DisjointTimeGrid::transmission(int /*ell*/, double /*time*/) const
{
    return 1.;
}

////////////////////////////////////////////////////////////////////

vector<int> DisjointTimeGrid::bins(double time) const
{
    // get the bin index, or -1 for "out of range"
    int ell = bin(time);

    // wrap the result in a list
    vector<int> result;
    if (ell >= 0) result.push_back(ell);
    return result;
}

////////////////////////////////////////////////////////////////////

int DisjointTimeGrid::bin(double time) const
{
    // get the index of the phantom time bin defined by the list of all K borders (where K=N+1 or K=N*2)
    //  0  => out of range on the left side
    //  K  => out of range on the right side
    size_t index = std::upper_bound(begin(_borderv), end(_borderv), time) - begin(_borderv);

    // map this index to the actual time bin index, or to -1 for "out of range"
    return _ellv[index];
}

////////////////////////////////////////////////////////////////////

Array DisjointTimeGrid::exttimev() const
{
    int n = _timev.size();
    Array extv(n + 2);

    extv[0] = _timeleftv[0];
    for (int ell = 0; ell != n; ++ell) extv[ell + 1] = _timev[ell];
    extv[n + 1] = _timerightv[n - 1];

    return extv;
}

////////////////////////////////////////////////////////////////////

Array DisjointTimeGrid::extdtimev() const
{
    int n = _dtimev.size();
    Array extv(n + 2);

    extv[0] = 0.;
    for (int ell = 0; ell != n; ++ell) extv[ell + 1] = _dtimev[ell];
    extv[n + 1] = 0.;

    return extv;
}

////////////////////////////////////////////////////////////////////
