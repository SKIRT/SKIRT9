/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef LOG_HPP
#define LOG_HPP

#include "SimulationItem.hpp"
#include <atomic>
#include <mutex>

////////////////////////////////////////////////////////////////////

/** Log is the abstract base class for a message logging mechanism. It offers convenience functions
    for logging messages at various levels (info, warning, success, error), adding a time-stamp
    along the way. All of these functions eventually call a single pure virtual function, which
    must be implemented in a subclass to actually output the message to a device such as the
    console or a file. The message-logging functions in this class are thread-safe. */
class Log : public SimulationItem
{
    //============= Construction - Setup - Destruction =============

protected:
    /** The default constructor is declared protected because this is an abstract class. */
    Log();

    /** The purpose of this function is setting the _procName attribute, using the rank of the process.
        This rank is obtained by using the find algorithm to search for an instance of a
        ProcessCommunicator subclass. If none is found, because for example the simulation or fitscheme
        hierarchy has not yet been created, the function returns and _procName retains its empty state.
        If the find algorithm does succeed, the find algorithm is applied again, this time to invoke
        the setup of the ProcessCommunicator object, if this is not yet done. Then, the rank is
        obtained from the ProcessCommunicator and passed as an argument to the setProcessName()
        function. */
    void setupSelfBefore() override;

private:
    /** This private function sets the process name used for logging based on the process rank
        passed as an argument to this function. This function makes sure that the same procedure is
        also applied for the linked Log instance, if present. */
    void setRank(int rank);

    //======== Setters & Getters for Discoverable Attributes =======

public:
    /** This enum includes a constant for each logging level, in increasing order of importance. */
    enum class Level { Info, Warning, Success, Error };

    /** Sets the lowest logging level that actually gets written to the output device. In other
        words, any messages logged at a lower level are ignored. The default value is Info, so that
        all messages are logged. */
    void setLowestLevel(Level level);

    /** Returns the lowest logging level that actually gets written to the output device. */
    Level lowestLevel() const;

    /** Sets the Log instance that is linked into this one. All messages received by this Log
        instance are also sent to the linked Log instance, regardless of the logging level set for
        this instance (i.e. each Log instance in the chain decides for itself which messages to
        display). The receiving Log instance assumes ownership of the linked Log instance. */
    void setLinkedLog(Log* log);

    /** Returns the Log instance that is linked into this one, so that it receives a copy of all
        messages. */
    Log* linkedLog() const;

    /** Sets or unsets the verbose mode for this Log instance. */
    void setVerbose(bool value);

    /** Returns whether the Log is set in verbose mode or not. */
    bool verbose() const;

    /** Sets or unsets the memory logging mode for this Log instance. */
    void setMemoryLogging(bool value);

    /** Returns whether memory usage is logged or not. */
    bool memoryLogging() const;

    //======================== Logging Functions =======================

public:
    /** Logs an informational message (i.e. at level Info). In multiprocessing mode, only the info
        messages of the root processes are actually logged. If verbose mode is enabled, however,
        all processes log and the process name is attached to the info message. This function is
        thread-safe. */
    void info(string message);

    /** Resets the interval timer for progress messages issued through the infoIfElapsed()
        function. The first argument specifies the total number of tasks to be performed, which
        will be used to calculate a completion percentage in each progress message. Use a value of
        zero to disable this feature. The second argument optionally specifies the minimum time
        interval between messages in seconds; the default is 5 seconds. Call this function once
        before issuing a sequence of potentially frequent progress messages through the
        infoIfElapsed() function. This function should be called only from a single thread (for a
        given Log instance). */
    void infoSetElapsed(size_t numTotal, int seconds = 5);

    /** Calls the info() function to log an informational message if a certain time interval
        (specified with the infoSetElapsed() function) has elapsed since the previous invocation of
        the infoIfElapsed() function (or of the infoSetElapsed() function). If the value of \em
        numTotal passed to infoSetElapsed() is nonzero, this function accumulates the values of \em
        numDone passed to the successive calls of infoIfElapsed() and uses this running subtotal to
        add a percentage of completion to the basic log message. This function is thread-safe. */
    void infoIfElapsed(string message, size_t numDone);

    /** Logs a warning message (i.e. at level Warning). Each warning message is logged,
        irrespective of which process invokes this function. The warning message is prefixed with
        the process name, if multiprocessing mode is used. If the second argument is true or
        omitted, the warning message is also stored for later retrieval (e.g., at the end of the
        simulation). This function is thread-safe. */
    void warning(string message, bool store = true);

    /** Logs an informational message (i.e. at level Success). In multiprocessing mode, only the
        success messages of the root processes are actually logged. If verbose mode is enabled,
        however, all processes log and the process name is attached to the success message. This
        function is thread-safe. */
    void success(string message);

    /** Logs an informational message (i.e. at level Error). Each error message is logged,
        irrespective of which process invokes this function. The error message is prefixed with the
        process name, if multiprocessing mode is used. If the second argument is true or omitted,
        the warning message is also stored for later retrieval (e.g., at the end of the
        simulation). This function is thread-safe. */
    void error(string message, bool store = true);

protected:
    /** This pure virtual function must be implemented in a subclass to actually output the
        specified message to a device such as the console or a file in a thread-safe way. The
        message does not yet contain a time stamp. The second argument specifies the logging level
        for the message (info, warning, success, error). The level is guaranteed to be at or above
        the current lowest level. */
    virtual void output(string message, Level level) = 0;

    /** This function returns a string identifying this process of the form "Pnnn",
        where nnn is the rank of the process. In singleprocessing mode, this string is empty. */
    string processName();

    //======================== Other Functions =======================

public:
    /** This function returns a list of the warning messages issued since this Log instance was
        constructed, limited to the first nine warnings, or the empty list if no warnings were
        issued. The messages are returned as they were passed to the warning() function, i.e.
        without time stamp, process name or other embellishment. This list can be used, for
        example, to repeat the warnings at the end of a simulation. */
    vector<string> warningsIssued();

    /** This function returns a list of the error messages issued since this Log instance was
        constructed, limited to the first nine errors, or the empty list if no errors were
        issued. The messages are returned as they were passed to the error() function, i.e.
        without time stamp, process name or other embellishment. This list can be used, for
        example, to repeat the errors at the end of a simulation. */
    vector<string> errorsIssued();

    //======================== Data Members ========================

private:
    Level _lowestLevel{Level::Info};
    Log* _link{nullptr};
    bool _verbose{false};
    bool _logmemory{false};
    string _procNameShort;
    string _procNameLong;
    std::atomic<uint64_t> _interval{0};  // counts in interval between messages
    std::atomic<uint64_t> _started{0};   // counts since epoch when timer started
    std::atomic<size_t> _numTotal{0};    // total number of tasks to be completed
    std::atomic<size_t> _numDone{0};     // number of tasks completed so far

    // list of warnings and errors issued, so they can be repeated at the end of the simulation
    vector<string> _warnings;
    vector<string> _errors;
    std::mutex _mutex;  // mutex to guard access to these message lists
};

////////////////////////////////////////////////////////////////////

#endif
