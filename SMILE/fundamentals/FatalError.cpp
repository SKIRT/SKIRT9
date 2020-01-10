/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "FatalError.hpp"
#include "System.hpp"
#include <sstream>

////////////////////////////////////////////////////////////////////

FatalError::FatalError(string message, const char* file, int line, const char* function)
{
    // split the message in lines if needed, and store the result
    std::stringstream ss(message);
    string msgline;
    while (getline(ss, msgline))
    {
        if (!msgline.empty()) _message.push_back(msgline);
    }

    // ensure that there is a first line "describing" the error
    if (_message.empty()) _message.push_back("Unknown error");

    // add information on the source code location
    _message.push_back("On line " + std::to_string(line) + " in file " + string(file));
    _message.push_back("In function " + string(function));

    // add a simple stack trace
    for (auto msgline : System::stacktrace()) _message.push_back(msgline);
}

////////////////////////////////////////////////////////////////////

vector<string> FatalError::message() const
{
    return _message;
}

////////////////////////////////////////////////////////////////////
