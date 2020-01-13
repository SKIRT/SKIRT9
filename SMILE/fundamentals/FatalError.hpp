/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef FATALERROR_HPP
#define FATALERROR_HPP

#include "Basics.hpp"

////////////////////////////////////////////////////////////////////

/** This class represents a fatal error which can be thrown as an exception from anywhere in the
    application, even during parallel execution. The top-level function in the application should
    catch the exception and alert the user, including the message held by the FatalError instance.
    An instance of FatalError should be thrown by value and caught by reference. */
class FatalError
{
public:
    /** Constructs a fatal error with the specified error message and information about the call
        site. The error message can be split over multiple lines by using newline characters
        (backslash n) as usual. Rather than directly using this constructor, one should use the
        FATALERROR macro defined in this header to automatically include the source code location
        information. For example:
        \code
        throw FATALERROR("Theta should be between 0 and pi.");
        \endcode
    */
    FatalError(string message, const char* file, int line, const char* function);

    /** Returns the multi-line error message for this fatal error as a list of strings. Each string
        contains a single line, i.e. the strings contain no newline characters. The actual error
        message as passed to the constructor is contained in the first line(s). Subsequent lines
        describe the location in the source code where the exception was thrown, and provide a
        basic call stack overview. */
    vector<string> message() const;

private:
    // the multi-line message, including source code location information and call stack overview.
    vector<string> _message;
};

////////////////////////////////////////////////////////////////////

/** Constructs a fatal error with the specified error message, and automatically includes the
    relevant information about the call site. */
#define FATALERROR(message) FatalError((message), __FILE__, __LINE__, __func__)

////////////////////////////////////////////////////////////////////

#endif
