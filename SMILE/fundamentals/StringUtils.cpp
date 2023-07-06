/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "StringUtils.hpp"
#include <locale>

////////////////////////////////////////////////////////////////////

bool StringUtils::startsWith(string text, string find)
{
    return text.length() >= find.length() && text.compare(0, find.length(), find) == 0;
}

////////////////////////////////////////////////////////////////////

bool StringUtils::endsWith(string text, string find)
{
    return text.length() >= find.length() && text.compare(text.length() - find.length(), find.length(), find) == 0;
}

////////////////////////////////////////////////////////////////////

bool StringUtils::contains(string text, string find)
{
    return text.find(find) != string::npos;
}

////////////////////////////////////////////////////////////////////

namespace
{
    // This function returns true if two given C strings match,
    // where the first string may contain the wildcard characters '*' and '?'.
    // The recursive implementation was copied from a web post by Vishal Chaudhary.
    bool match(const char* first, const char* second)
    {
        // If we reach at the end of both strings, we are done
        if (*first == '\0' && *second == '\0') return true;

        // Make sure that the characters after '*' are present
        // in second string. This function assumes that the first
        // string will not contain two consecutive '*'
        if (*first == '*' && *(first + 1) != '\0' && *second == '\0') return false;

        // If the first string contains '?', or current characters
        // of both strings match, test the remaining segment
        if (*first == '?' || *first == *second) return match(first + 1, second + 1);

        // If there is *, then there are two possibilities
        // a) We consider current character of second string
        // b) We ignore current character of second string.
        if (*first == '*') return match(first + 1, second) || match(first, second + 1);
        return false;
    }
}

////////////////////////////////////////////////////////////////////

bool StringUtils::matches(string text, string pattern)
{
    return match(pattern.c_str(), text.c_str());
}

////////////////////////////////////////////////////////////////////

bool StringUtils::contains(const vector<string>& list, string find)
{
    return std::find(list.cbegin(), list.cend(), find) != list.cend();
}

////////////////////////////////////////////////////////////////////

int StringUtils::indexOf(const vector<string>& list, string find)
{
    auto it = std::find(list.cbegin(), list.cend(), find);
    if (it != list.cend())
        return static_cast<int>(it - list.cbegin());
    else
        return -1;
}

////////////////////////////////////////////////////////////////////

string StringUtils::replace(string text, string find, string replace)
{
    for (string::size_type i = 0; (i = text.find(find, i)) != string::npos; i += replace.length())
    {
        text.replace(i, find.length(), replace);
    }
    return text;
}

////////////////////////////////////////////////////////////////////

namespace
{
    bool isWhiteSpace(char c) { return c == ' ' || c == '\t' || c == 0x0D || c == 0x0A; }
}

string StringUtils::squeeze(string text)
{
    auto destination = text.begin();
    bool haveSpace = true;

    for (auto source = text.cbegin(); source != text.cend(); ++source)
    {
        if (isWhiteSpace(*source))
        {
            if (!haveSpace)
            {
                *destination++ = ' ';
                haveSpace = true;
            }
        }
        else
        {
            *destination++ = *source;
            haveSpace = false;
        }
    }
    if (haveSpace && destination != text.cbegin()) destination--;

    text.erase(destination, text.end());
    return text;
}

////////////////////////////////////////////////////////////////////

string StringUtils::padLeft(string text, size_t length, char pad)
{
    if (text.length() < length) text.insert(0, length - text.length(), pad);
    return text;
}

////////////////////////////////////////////////////////////////////

string StringUtils::padRight(string text, size_t length, char pad)
{
    if (text.length() < length) text.append(length - text.length(), pad);
    return text;
}

////////////////////////////////////////////////////////////////////

string StringUtils::toLower(string text)
{
    // use the C++ function, which works with any character type (as opposed to the less portable C function)
    // use the global C++ locale, which is set to the standard "C" locale in System::initialize()
    std::locale loc;
    for (auto& c : text)  // writable reference to character in string
    {
        c = std::tolower(c, loc);
    }
    return text;
}

////////////////////////////////////////////////////////////////////

string StringUtils::toUpper(string text)
{
    // use the C++ function, which works with any character type (as opposed to the less portable C function)
    // use the global C++ locale, which is set to the standard "C" locale in System::initialize()
    std::locale loc;
    for (auto& c : text)  // writable reference to character in string
    {
        c = std::toupper(c, loc);
    }
    return text;
}

////////////////////////////////////////////////////////////////////

string StringUtils::toUpperFirst(string text)
{
    if (!text.empty())
    {
        std::locale loc;
        text[0] = std::toupper(text[0], loc);
    }
    return text;
}

////////////////////////////////////////////////////////////////////

vector<string> StringUtils::split(string text, string separator)
{
    vector<string> result;

    string::size_type start = 0;
    string::size_type end = text.find(separator);
    while (end != std::string::npos)
    {
        result.push_back(text.substr(start, end - start));
        start = end + separator.length();
        end = text.find(separator, start);
    }
    result.push_back(text.substr(start, end));
    return result;
}

////////////////////////////////////////////////////////////////////

string StringUtils::join(const vector<string>& list, string separator)
{
    string result;
    for (const string& segment : list)
    {
        if (!result.empty() && !segment.empty()) result += separator;
        result += segment;
    }
    return result;
}

////////////////////////////////////////////////////////////////////

string StringUtils::joinPaths(string segment1, string segment2)
{
    if (segment1.empty()) return segment2;
    if (segment2.empty()) return segment1;

    bool slash1 = endsWith(segment1, "/") || endsWith(segment1, "\\");
    bool slash2 = startsWith(segment2, "/") || startsWith(segment2, "\\");

    if (slash1 & slash2) return segment1 + segment2.substr(1);
    if (slash1 | slash2) return segment1 + segment2;
    return segment1 + "/" + segment2;
}

////////////////////////////////////////////////////////////////////

bool StringUtils::isAbsolutePath(string filepath)
{
    return startsWith(filepath, "/") || startsWith(filepath, "\\") || contains(filepath, ":");
}

////////////////////////////////////////////////////////////////////

string StringUtils::filename(string filepath)
{
    // remove trailing slash, if any
    if (filepath.length() > 1 && (filepath.back() == '/' || filepath.back() == '\\')) filepath.pop_back();

    // locate rightmost slash, if any, and return everything after that position
    auto index = filepath.find_last_of("/\\");
    if (index != string::npos) return filepath.substr(index + 1);

    // if there are no slashes, return the complete path
    return filepath;
}

////////////////////////////////////////////////////////////////////

string StringUtils::filenameBase(string filepath)
{
    // get the filename part, including extension
    filepath = filename(filepath);

    // locate rightmost period and remove everything from that position onwards
    auto index = filepath.rfind('.');
    if (index != string::npos) filepath.erase(index);

    return filepath;
}

////////////////////////////////////////////////////////////////////

string StringUtils::dirPath(string filepath)
{
    if (filepath.length() > 1 && filepath.find_first_of("/\\") != string::npos)
    {
        // remove trailing slash
        if (filepath.back() == '/' || filepath.back() == '\\') filepath.pop_back();

        // locate rightmost slash and remove everything after that position
        auto index = filepath.find_last_of("/\\");
        if (index != string::npos)
        {
            filepath.erase(index + 1);
            return filepath;
        }
    }
    return string();
}

////////////////////////////////////////////////////////////////////

string StringUtils::addExtension(string filename, string extension)
{
    if (!filename.empty() && !endsWith(StringUtils::toLower(filename), "." + toLower(extension)))
    {
        filename += "." + extension;
    }
    return filename;
}

////////////////////////////////////////////////////////////////////

bool StringUtils::isValidBool(string value)
{
    value = toLower(squeeze(value));
    if (value == "true" || value == "t" || value == "yes" || value == "y" || value == "1" || value == "false"
        || value == "f" || value == "no" || value == "n" || value == "0")
        return true;
    return false;
}

////////////////////////////////////////////////////////////////////

bool StringUtils::toBool(string value)
{
    value = toLower(squeeze(value));
    if (value == "true" || value == "t" || value == "yes" || value == "y" || value == "1") return true;
    return false;
}

////////////////////////////////////////////////////////////////////

string StringUtils::toString(bool value)
{
    return value ? "true" : "false";
}

////////////////////////////////////////////////////////////////////

bool StringUtils::isValidInt(string value)
{
    value = squeeze(value);
    if (value.empty()) return false;

    const char* begptr = value.c_str();
    char* endptr;
    auto result = std::strtol(begptr, &endptr, 10);
    return static_cast<string::size_type>(endptr - begptr) == value.length() && result >= INT_MIN && result <= INT_MAX;
}

////////////////////////////////////////////////////////////////////

int StringUtils::toInt(string value)
{
    return isValidInt(value) ? std::stoi(value) : 0;
}

////////////////////////////////////////////////////////////////////

string StringUtils::toString(int value)
{
    return std::to_string(value);
}

////////////////////////////////////////////////////////////////////

bool StringUtils::isValidDouble(string value)
{
    bool ok;
    toDouble(value, &ok);
    return ok;
}

////////////////////////////////////////////////////////////////////

double StringUtils::toDouble(string value, bool* ok)
{
    value = squeeze(value);
    if (!value.empty())
    {
        const char* begptr = value.c_str();
        char* endptr;
        auto result = std::strtod(begptr, &endptr);
        if (static_cast<string::size_type>(endptr - begptr) == value.length() && abs(result) < HUGE_VAL)
        {
            if (ok) *ok = true;
            return result;
        }
    }
    if (ok) *ok = false;
    return 0.;
}

////////////////////////////////////////////////////////////////////

string StringUtils::toString(double value)
{
    // use a decent representation for not-a-number and infinity
    if (std::isnan(value)) return "∅";
    if (std::isinf(value)) return value < 0 ? "-∞" : "∞";

    // start with a regular semi-smart conversion
    char buf[20];
    snprintf(buf, sizeof(buf), "%1.10g", value);
    string result(buf);

    // remove leading zeroes and the + sign in the exponent
    result = replace(result, "e-0", "e-");
    result = replace(result, "e+0", "e");
    result = replace(result, "e+", "e");

    // replace 4 or more trailing zeroes by exponent
    string::size_type zeroes = 0;
    for (auto it = result.crbegin(); it != result.crend(); ++it, ++zeroes)
        if (*it != '0') break;
    if (zeroes > 3)
    {
        result.erase(result.length() - zeroes);
        result += "e" + std::to_string(zeroes);
    }
    return result;
}

////////////////////////////////////////////////////////////////////

string StringUtils::toString(double value, char format, int precision, int width, char pad)
{
    // force 'd' and unknown formats to fixed point with zero digits after decimal point
    if (format != 'f' && format != 'e' && format != 'g')
    {
        format = 'f';
        precision = 0;
    }

    // avoid very large values in fixed point format (because they would overrun our output buffer)
    if (format == 'f' && value >= 1e20)
    {
        format = 'e';
        precision = 18;
    }

    // force precision within range
    precision = min(max(precision, 0), 18);

    // perform the conversion
    char formatString[] = {'%', '1', '.', '*', format, 0};
    char result[30];
    snprintf(result, sizeof(result), formatString, precision, value);

    // pad if needed
    return padLeft(result, static_cast<size_t>(max(width, 1)), pad);
}

////////////////////////////////////////////////////////////////////

string StringUtils::toMemSizeString(size_t value)
{
    static const char prefixes[] = {'K', 'M', 'G', 'T'};
    size_t factor = 1024;
    size_t index = 0;
    for (; index < 3; ++index)
    {
        if (value < factor * 999) break;
        factor *= 1024;
    }
    return toString(static_cast<double>(value) / static_cast<double>(factor), 'g', 3) + ' ' + prefixes[index] + 'B';
}

////////////////////////////////////////////////////////////////////
