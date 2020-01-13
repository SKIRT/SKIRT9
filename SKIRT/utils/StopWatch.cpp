/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "StopWatch.hpp"
#include "FatalError.hpp"
#include "StringUtils.hpp"
#include <chrono>

////////////////////////////////////////////////////////////////////

namespace
{
    const int N = 5;  // the number of timers exposed to the user
    int _level = -1;  // the current level == index in timer arrays

    // these arrays have an element for each timer
    uint64_t _count[N] = {0};  // total number of calls to stop()
    uint64_t _total[N] = {0};  // accumulated elapsed time between start/stop (in system ticks)
    uint64_t _start[N] = {0};  // absolute time at most recent start (in system ticks)

    // returns the current time as the number of system ticks gone by since some fixed reference time point
    uint64_t now()
    {
        using namespace std::chrono;
        return high_resolution_clock::now().time_since_epoch().count();
    }

    // returns the conversion factor from system ticks to seconds
    constexpr double tick()
    {
        using namespace std::chrono;
        return static_cast<double>(high_resolution_clock::period::num)
               / static_cast<double>(high_resolution_clock::period::den);
    }
}

////////////////////////////////////////////////////////////////////

StopWatch::StopWatch()
{
    start();
}

////////////////////////////////////////////////////////////////////

StopWatch::~StopWatch()
{
    stop();
}

////////////////////////////////////////////////////////////////////

void StopWatch::start()
{
    _level++;
    if (_level >= N) throw FATALERROR("Timer nesting overflow");

    _start[_level] = now();
}

////////////////////////////////////////////////////////////////////

void StopWatch::stop()
{
    if (_level < 0) throw FATALERROR("Timer nesting underflow");

    uint64_t stop = now();
    _count[_level]++;
    _total[_level] += stop - _start[_level];
    _start[_level] = stop;

    _level--;
}

////////////////////////////////////////////////////////////////////

vector<string> StopWatch::report()
{
    if (_level != -1) throw FATALERROR("Timer nesting inbalance");
    vector<string> result;

    // only produce a non-empty result if at least one timer was actually used
    if (_count[0])
    {
        // get the conversion factor from system ticks to seconds
        const double conv = tick();

        // add a line per timer that was actually used,
        // and calculate the total and maximum number of start/stop sequences
        uint64_t totalcount = 0;
        uint64_t maxcount = 0;
        double total0 = conv * _total[0];
        for (int i = 0; i < N; i++)
        {
            if (_count[i])
            {
                double total = conv * _total[i];
                result.push_back("Timer " + StringUtils::toString(i + 1) + ":"
                                 + StringUtils::toString(total, 'f', 3, 10) + " s  "
                                 + StringUtils::toString(100 * total / total0, 'f', 1, 5) + " %");
                totalcount += _count[i];
                maxcount = max(maxcount, _count[i]);
            }
        }

        // estimate the start/stop time using the first timer
        const int K = 5000;  // number of start/stop sequences in the test
        _total[0] = 0;
        for (int k = 0; k < K; k++)
        {
            start();
            stop();
        }
        double startstop = conv * _total[0] / K;

        // add a line with error information
        // estimate the total error: start/stop time accumulates across timers, measurement resolution doesn't
        double error = startstop * totalcount + conv * maxcount;
        result.push_back("Error ±:" + StringUtils::toString(error, 'f', 3, 10) + " s  "
                         + StringUtils::toString(100 * error / total0, 'f', 1, 5) + " %");
    }
    return result;
}

////////////////////////////////////////////////////////////////////
