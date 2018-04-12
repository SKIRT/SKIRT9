/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef CONSOLELOG_HPP
#define CONSOLELOG_HPP

#include "Log.hpp"

////////////////////////////////////////////////////////////////////

/** ConsoleLog inherits from Log and implements logging to the standard console. All ConsoleLog
    instances share the same underlying console device. It is safe to mix multiple instances. The
    output() function in this class is thread-safe. */
class ConsoleLog : public Log
{
    //============= Construction - Setup - Destruction =============

public:
    /** This constructor creates a standalone console log object that is not hooked into a
        simulation item hierarchy; the caller is responsible for its destruction. The setup()
        function is \em not called by this constructor. The standalone object can be used without
        further setup, however the generated messages will contain no multiprocessing information
        (such as the process rank) until after setup() has been called. */
    ConsoleLog();

    /** This constructor creates a console log object that is hooked up as a child to the specified
        parent in the simulation hierarchy, so that it will automatically be deleted. The setup()
        function is \em not called by this constructor. */
    explicit ConsoleLog(SimulationItem* parent);

    //======================== Other Functions =======================

protected:
    /** This function outputs a message to the console, colored and annotated according to the
        specified logging level. It overrides the pure virtual function in the base class. This
        function is thread-safe. */
    void output(string message, Level level) override;
};

////////////////////////////////////////////////////////////////////

#endif
