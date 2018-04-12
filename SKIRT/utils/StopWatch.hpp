/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef STOPWATCH_HPP
#define STOPWATCH_HPP

#include "Basics.hpp"

////////////////////////////////////////////////////////////////////

/**
This class enables timing of code snippets by summing the elapsed time for each pass through
the snippet. It uses the standard C++ high-resolution clock for maximum precision. The accuracy
of the results is platform-specific, as it depends on the resolution of the underlying time
measurement and by the time needed to start and stop the watch.

The functions in this class are \em not thread-safe and \em not re-entrant. The results are
meaningful only if the application uses just a single execution thread and no other processes
consume time on the computer. Even then interpreting timings is tricky. For example, some
modern desktop systems dynamically change the CPU clock speed depending on the processor's chip
temperature, and the placement of data in memory blocks relative to the processor core in use
may substantially affect performance.

This class offers 5 global (i.e. application-wide) timers numbered from 1 to 5. The timers are
nested, i.e. timer n+1 can only be started when timer n is already running. Vice versa, timer n
can only be stopped when timer n+1 is no longer running. At the start of the application all
timers are initialized to zero and the nesting level is set to zero as well. A call to the
static start() function increments the nesting level and starts the timer corresponding to the
new level. A call to the static stop() function stops the timer corresponding to the current
level and decrements the nesting level. If the nesting level goes out of range, a fatal error
is thrown. To access the timing results, the static report() function returns a list of strings
reporting information on all nonzero timers in a human-readable format.

Rather than explicitly invoking the static start() and stop() functions, an instance of the
StopWatch class can be used to ensure correct nesting of the timings: the constructor simply
calls start() and the destructor calls stop(). For example, to time the execution of a complete
scope, regardless of the scope's exit point, construct a StopWatch instance at the start of the
scope:

\code
{
    StopWatch watch;
    ...
    if (...) return 1;
    else if (...) return 2;
    ...
}
\endcode

*/
class StopWatch
{
public:
    /** The constructor calls the static start() function to start the timer at the next nesting
        level. */
    StopWatch();

    /** The destructor calls the static stop() function to stop the timer at the current nesting
        level. */
    ~StopWatch();

    /** This function increments the nesting level and starts the timer corresponding to the new
        level. If the nesting level goes out of range, a fatal error is thrown. */
    static void start();

    /** This function stops the timer corresponding to the current level and decrements the nesting
        level. If the nesting level goes out of range, a fatal error is thrown.*/
    static void stop();

    /** This function returns a list of strings reporting information on all nonzero timers in a
        human-readable format. When this function is called, none of the timers may be running.
        If all timers are still zero, the function returns an empty list. */
    static vector<string> report();
};

////////////////////////////////////////////////////////////////////

#endif
