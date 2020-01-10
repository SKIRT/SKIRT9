/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "Console.hpp"
#include "FatalError.hpp"
#include "StringUtils.hpp"
#include "System.hpp"

////////////////////////////////////////////////////////////////////

void Console::info(string message)
{
    System::log(message, System::LogLevel::Info);
}

////////////////////////////////////////////////////////////////////

void Console::warning(string message)
{
    System::log(message, System::LogLevel::Warning);
}

////////////////////////////////////////////////////////////////////

void Console::success(string message)
{
    System::log(message, System::LogLevel::Success);
}

////////////////////////////////////////////////////////////////////

void Console::error(string message)
{
    System::log(message, System::LogLevel::Error);
}

////////////////////////////////////////////////////////////////////

bool Console::promptForBool(string message, bool hasDef, bool def)
{
    // add hints on possible and default values to the message
    message += " [yes/no]";
    if (hasDef) message += def ? " (yes)" : " (no)";

    while (true)
    {
        // get the input string
        string input = System::prompt(message);

        // if provided, use the default value instead of an empty string
        if (input.empty() && hasDef) return def;

        // if successful conversion, return result
        if (StringUtils::isValidBool(input)) return StringUtils::toBool(input);

        // reject conversion errors or empty strings without default
        error("Enter 'yes' or 'no'");
    }
}

////////////////////////////////////////////////////////////////////

int Console::promptForInt(string message, int min, int max, bool hasDef, int def)
{
    // verify that default is in range
    if (hasDef && (def < min || def > max)) throw FATALERROR("Default value out of range");

    // add hints on min, max and default values to the message
    message += " [" + StringUtils::toString(min) + "," + StringUtils::toString(max) + "]";
    if (hasDef) message += " (" + StringUtils::toString(def) + ")";

    while (true)
    {
        // get the input string and attempt conversion
        string input = System::prompt(message);
        int result = StringUtils::toInt(input);

        // if provided, use the default value instead of an empty string
        if (input.empty() && hasDef) return def;

        // reject conversion errors or empty strings without default
        if (!StringUtils::isValidInt(input))
        {
            error("Enter a valid integer number");
        }

        // reject out-of-range values
        else if (result < min)
        {
            error("Enter a number larger than or equal to " + StringUtils::toString(min));
        }
        else if (result > max)
        {
            error("Enter a number smaller than or equal to " + StringUtils::toString(max));
        }
        // return successful result
        else
        {
            return result;
        }
    }
}

////////////////////////////////////////////////////////////////////

string Console::promptForString(string message, bool hasDef, string def)
{
    // add hint on default value to the message
    if (hasDef) message += " (" + def + ")";

    while (true)
    {
        // get the input string
        string input = System::prompt(message);

        // accept non-empty string
        if (!input.empty()) return input;

        // if provided, use the default value instead of an empty string
        if (hasDef) return def;

        // reject empty strings
        error("Enter a nonempty string");
    }
}

////////////////////////////////////////////////////////////////////

int Console::promptForChoice(string message, const vector<string>& choices, bool hasDef, int defIndex,
                             bool allowNoChoice, string noChoiceMessage)
{
    if (choices.empty() && !allowNoChoice) throw FATALERROR("There are no choices to prompt for");
    if (defIndex < 0) hasDef = false;

    info("Possible choices for " + message + ":");
    for (size_t index = 0; index < choices.size(); index++)
    {
        string choice = StringUtils::toUpperFirst(choices[index]);
        info(string(index < 9 ? "   " : "  ") + std::to_string(index + 1) + ". " + choice);
    }

    if (choices.size() == 1 && !allowNoChoice)
    {
        info("Automatically selected the only choice: 1");
        return 0;
    }
    if (choices.empty() && allowNoChoice)
    {
        info("Automatically selected the only choice: 0");
        return -1;
    }
    return promptForInt("Enter one of these numbers" + (allowNoChoice ? " " + noChoiceMessage : ""),
                        (allowNoChoice ? 0 : 1), static_cast<int>(choices.size()), true,
                        (hasDef ? defIndex + 1 : (allowNoChoice ? 0 : 1)))
           - 1;
}

////////////////////////////////////////////////////////////////////
