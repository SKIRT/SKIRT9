/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef COMMANDLINEARGUMENTS_HPP
#define COMMANDLINEARGUMENTS_HPP

#include "Basics.hpp"
#include <unordered_map>

////////////////////////////////////////////////////////////////////

/** This class serves as a basic command line parser. It supports any number of filepath arguments
    and options (with or without a value). Individual filepaths, options, and option values must be
    separated from each other by whitespace. An option always starts with a dash; a filepath and an
    option value don't. Filepaths and options can occur in any (mixed) order, except that an option
    value should immediately follow the corresponding option. The filepath arguments can be
    retrieved in order of occurrence. The presence of an option (and if applicable its value) can
    be retrieved by option name. If an option occurs multiple times, the last value is retained.
    Other than that, the order of options is not preserved. */
class CommandLineArguments final
{
public:
    /** This constructor parses the specified list of command line arguments, usually retrieved
        through the System::arguments() function, according to the specified option list, and
        stores the result internally. The option list is specified as a string containing the
        allowed options separated by whitespace, in arbitrary order. Options taking a value should
        be followed by an asterisk (without whitespace); the asterisk is removed from the option
        name. For example, the option list "-t* -o* -b -opt" means that there are four allowed
        options: -t and -o take a value, while -b and -opt don't. Filepath arguments are not
        specified as they are always allowed before, after, or mixed in with the options. */
    CommandLineArguments(const vector<string>&, string options);

    /** Returns true if the command line is valid, i.e. if it contains only allowed options and a
        value is provided for options that take a value. */
    bool isValid() const;

    /** Returns true if there is at least one option (with or without value), or false if not. If
        the command line is invalid, this function always returns false. */
    bool hasOptions() const;

    /** Returns true if the specified option is present (with or without a value), or false if not.
        If the command line is invalid, this function always returns false. */
    bool isPresent(string option) const;

    /** Returns the value of the specified option, or the empty string if the option is not present
        or if the option does not take a value. If the command line is invalid, this function
        always returns an empty string. */
    string value(string option) const;

    /** Returns the value of the specified option converted to an integer, or -1 if the option is
        not present or if the value can't be converted to an integer. If the command line is
        invalid, this function always returns -1. Since an option value string can't start with a
        dash, it's impossible to represent negative integers, thus an error return value of -1 is
        unambiguous. */
    int intValue(string option) const;

    /** Returns the value of the specified option converted to a double value, or -1 if the option
        is not present or if the value can't be converted to a double. If the command line is
        invalid, this function always returns -1. Since an option value string can't start with a
        dash, it's impossible to represent negative doubles, thus an error return value of -1 is
        unambiguous. */
    double doubleValue(string option) const;

    /** Returns true if there is at least one filepath argument, or false if not. If the command
        line is invalid, this function always returns false. */
    bool hasFilepaths() const;

    /** Returns a list of the filepath arguments, in order of occurrence on the command line. If
        there are no filepath arguments, or if the command line is invalid, this function returns
        an empty list. */
    vector<string> filepaths() const;

private:
    // data members
    bool _valid;
    std::unordered_map<string, string> _options;
    vector<string> _filepaths;
};

////////////////////////////////////////////////////////////////////

#endif
