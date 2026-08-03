/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "TimeLogger.hpp"
#include "Log.hpp"
#include "StringUtils.hpp"
#include <exception>

////////////////////////////////////////////////////////////////////

TimeLogger::TimeLogger(Log* log, string scope) : _log(log), _scope(scope), _started(std::chrono::steady_clock::now())
{
    if (log) log->info("Starting " + scope + "...");
}

////////////////////////////////////////////////////////////////////

namespace
{
    // local constants
    const int64_t msecsInSecond = 1000;
    const int64_t msecsInMinute = 60 * msecsInSecond;
    const int64_t msecsInHour = 60 * msecsInMinute;
    const int64_t msecsInDay = 24 * msecsInHour;
}

////////////////////////////////////////////////////////////////////

TimeLogger::~TimeLogger()
{
    using namespace std::chrono;

    // If no Log instance was passed, we don't have to calculate the elapsed time
    if (!_log) return;

    // if this destructor is executed while an exception is unwinding the stack,
    // then we shouldn't report a success message!
    if (std::uncaught_exception()) return;

    // get the elapsed time in milliseconds
    int64_t msecs = duration_cast<milliseconds>(steady_clock::now() - _started).count();
    bool moreThanMinute = msecs >= msecsInMinute;

    // always include the elapsed time in seconds
    string elapsed = StringUtils::toString(static_cast<double>(msecs) / static_cast<double>(msecsInSecond), 'f',
                                           moreThanMinute ? 0 : 1)
                     + " s";

    // if the elapsed time is over a minute, also include a "0d 0h 0m 0s" format
    if (moreThanMinute)
    {
        int64_t days = msecs / msecsInDay;
        msecs -= days * msecsInDay;
        int64_t hours = msecs / msecsInHour;
        msecs -= hours * msecsInHour;
        int64_t minutes = msecs / msecsInMinute;
        msecs -= minutes * msecsInMinute;
        int64_t seconds = (msecs + msecsInSecond / 2) / msecsInSecond;

        elapsed += " (";
        if (days) elapsed += std::to_string(days) + "d ";
        if (days || hours) elapsed += std::to_string(hours) + "h ";
        elapsed += std::to_string(minutes) + "m ";
        elapsed += std::to_string(seconds) + "s)";
    }

    _log->success("Finished " + _scope + " in " + elapsed + ".");
}

////////////////////////////////////////////////////////////////////
