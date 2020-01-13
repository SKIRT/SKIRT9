/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "CommandLineArguments.hpp"
#include "StringUtils.hpp"

////////////////////////////////////////////////////////////////////

CommandLineArguments::CommandLineArguments(const vector<string>& cmdlineargs, string options)
{
    // parse the option list into a dictionary with a value of true if the option takes a value
    std::unordered_map<string, bool> takesValue;
    for (string option : StringUtils::split(StringUtils::squeeze(options), " "))
    {
        // ignore options that don't start with a dash
        if (StringUtils::startsWith(option, "-"))
        {
            if (StringUtils::endsWith(option, "*"))
                takesValue[option.substr(0, option.length() - 1)] = true;
            else
                takesValue[option] = false;
        }
    }

    // parse the command line arguments and verify allowed options and values, aborting when invalid
    _valid = false;
    for (size_t index = 0; index < cmdlineargs.size(); ++index)
    {
        string arg = cmdlineargs[index];

        // process an option and, if applicable, its value
        if (StringUtils::startsWith(arg, "-"))
        {
            if (!takesValue.count(arg)) return;
            if (takesValue[arg])
            {
                if (index + 1 >= cmdlineargs.size()) return;
                string value = cmdlineargs[++index];
                if (StringUtils::startsWith(value, "-")) return;
                _options[arg] = value;
            }
            else
                _options[arg] = "";
        }

        // process a filepath
        else
        {
            _filepaths.push_back(arg);
        }
    }
    _valid = true;
}

////////////////////////////////////////////////////////////////////

bool CommandLineArguments::isValid() const
{
    return _valid;
}

////////////////////////////////////////////////////////////////////

bool CommandLineArguments::hasOptions() const
{
    return _valid && !_options.empty();
}

////////////////////////////////////////////////////////////////////

bool CommandLineArguments::isPresent(string option) const
{
    return _valid && _options.count(option);
}

////////////////////////////////////////////////////////////////////

string CommandLineArguments::value(string option) const
{
    return _valid && _options.count(option) ? _options.at(option) : "";
}

////////////////////////////////////////////////////////////////////

int CommandLineArguments::intValue(string option) const
{
    string stringvalue = value(option);
    return StringUtils::isValidInt(stringvalue) ? StringUtils::toInt(stringvalue) : -1;
}

////////////////////////////////////////////////////////////////////

double CommandLineArguments::doubleValue(string option) const
{
    string stringvalue = value(option);
    return StringUtils::isValidDouble(stringvalue) ? StringUtils::toDouble(stringvalue) : -1;
}

////////////////////////////////////////////////////////////////////

bool CommandLineArguments::hasFilepaths() const
{
    return _valid && !_filepaths.empty();
}

////////////////////////////////////////////////////////////////////

vector<string> CommandLineArguments::filepaths() const
{
    return _valid ? _filepaths : vector<string>();
}

////////////////////////////////////////////////////////////////////
