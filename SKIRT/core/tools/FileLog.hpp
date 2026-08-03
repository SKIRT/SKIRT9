/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef FILELOG_HPP
#define FILELOG_HPP

#include "Log.hpp"
#include <fstream>
#include <mutex>

////////////////////////////////////////////////////////////////////

/** FileLog inherits from Log and implements thread-safe logging to a file. The file has a fixed
    name <tt>prefix_log.txt</tt> and is placed in the output filepath provided by the FilePaths
    instance attached to the simulation hierarchy at setup time. The log text is written in UTF-8
    encoding. The output() function in this class is thread-safe. */
class FileLog : public Log
{
    //============= Construction - Setup - Destruction =============

public:
    /** This constructor creates a file log object that is not hooked into a simulation item
        hierarchy, and has not been setup. The caller is responsible for calling the setup()
        function and for destroying the object. Alternatively, the object can be linked into a
        simulation item hierarchy by passing it to the setLinkedLog() function of another log
        object. */
    FileLog();

    /** This constructor creates a file log object that is hooked up as a child to the specified
        parent in the simulation hierarchy, so that it will automatically be deleted. The setup()
        function is \em not called by this constructor. */
    explicit FileLog(SimulationItem* parent);

protected:
    /** This function constructs the filename and opens the log file, overwriting any existing file
        with the same name. */
    void setupSelfBefore() override;

private:
    /** This function provides the implementation of opening the file, called by setupSelfBefore().
        With multiple processes and when not in verbose mode, this function can be called later on
        when a warning or error is encountered on one of the processes. */
    void open();

    //======================== Other Functions =======================

protected:
    /** This function outputs a message to the file. It overrides the pure virtual function in the
        base class. This function is thread-safe. */
    void output(string message, Level level) override;

    //======================== Data Members ========================

private:
    std::mutex _mutex;  // mutex to guard the input/output operations
    std::ofstream _out;
};

////////////////////////////////////////////////////////////////////

#endif
