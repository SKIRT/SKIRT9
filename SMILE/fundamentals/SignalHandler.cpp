/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "SignalHandler.hpp"
#include "FatalError.hpp"
#include <signal.h>

//////////////////////////////////////////////////////////////////////

// this function can't be in a local namespace (it compiles but does not work)
static void my__signal__handler(int sig)
{
    // to avoid indefinite loops: if we arrive here a second time, simply quit the application
    static bool haveDoneThis = false;
    if (haveDoneThis) exit(2);
    haveDoneThis = true;

    // throw an error with an appropriate message
    switch (sig)
    {
        case SIGABRT: throw FATALERROR("SIGNAL: Abort");
        case SIGFPE: throw FATALERROR("SIGNAL: Floating Point Exception");
        case SIGILL: throw FATALERROR("SIGNAL: Illegal Instruction");
        case SIGSEGV: throw FATALERROR("SIGNAL: Segmentation Violation");
        default: throw FATALERROR("SIGNAL: Unknown");
    }
}

////////////////////////////////////////////////////////////////////

void SignalHandler::InstallSignalHandlers()
{
    signal(SIGABRT, my__signal__handler);
    signal(SIGFPE, my__signal__handler);
    signal(SIGILL, my__signal__handler);
    signal(SIGSEGV, my__signal__handler);
}

////////////////////////////////////////////////////////////////////
