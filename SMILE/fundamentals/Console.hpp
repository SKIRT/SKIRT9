/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef CONSOLE_HPP
#define CONSOLE_HPP

#include "Basics.hpp"

////////////////////////////////////////////////////////////////////

/** The Console class offers facilities for interacting with the user through the console on a
    higher level than the basic functions in the System class. There is a seperate function to log
    a message at each log level, and there are functions to prompt for various data types,
    including a choice from a list, with proper validation. */
class Console final
{
public:
    // ================== Logging ==================

    /** Logs an informational message (i.e. at level Info). */
    static void info(string message);

    /** Logs a warning message (i.e. at level Warning). */
    static void warning(string message);

    /** Logs an informational message (i.e. at level Success). */
    static void success(string message);

    /** Logs an informational message (i.e. at level Error). */
    static void error(string message);

    // ================== Prompting ==================

    /** Prompts the console for a yes/no reply within the specified default value, and returns the
        user's response. */
    static bool promptForBool(string message, bool hasDef, bool def);

    /** Prompts the console for an integer number within the specified range and with the specified
        default value, and returns the user's response. */
    static int promptForInt(string message, int min, int max, bool hasDef, int def);

    /** Prompts the console for a non-empty string with the specified default value, and returns
        the user's response. */
    static string promptForString(string message, bool hasDef, string def);

    /** Prompts the console for a choice from the specified list, with the specified default, and
        returns the user's response. The function returns a zero-based index into the \em choices
        list. If \em allowNoChoice is true, the function returns -1 to indicate that no choice was
        made. */
    static int promptForChoice(string message, const vector<string>& choices, bool hasDef = false, int defIndex = -1,
                               bool allowNoChoice = false, string noChoiceMessage = string());
};

////////////////////////////////////////////////////////////////////

#endif
