/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef SIGNALHANDLER_HPP
#define SIGNALHANDLER_HPP

#include "Basics.hpp"

////////////////////////////////////////////////////////////////////

/** The SignalHandler class offers facilities for installing handlers for system errors. */
class SignalHandler final
{
public:
    /** This function installs C signal handlers for the error signals that may occur during
        execution, including such horrors as segmentation violation and bus error. The handlers
        throw a FatalError with an appropriate message, which provides a stack trace (assuming that
        the state of the program is not too messed up for memory allocation and such to work). */
    static void InstallSignalHandlers();
};

////////////////////////////////////////////////////////////////////

#endif
