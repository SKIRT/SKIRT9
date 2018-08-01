/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef STRINGUTILS_HPP
#define STRINGUTILS_HPP

#include "Basics.hpp"

////////////////////////////////////////////////////////////////////

/** The StringUtils class offers various facilities for working with std::string objects. */
class StringUtils final
{
public:
    // ================== Queries ==================

    /** Returns true if the text string starts with the find string, or if the find string is
        empty. Otherwise returns false. */
    static bool startsWith(string text, string find);

    /** Returns true if the text string ends with the find string, or if the find string is
        empty. Otherwise returns false. */
    static bool endsWith(string text, string find);

    /** Returns true if the text string contains the find string. */
    static bool contains(string text, string find);

    /** Returns true if the text string matches the pattern string, which may include * and ? characters. */
    static bool matches(string text, string pattern);

    /** Returns true if the string list contains the find string. */
    static bool contains(const vector<string>& list, string find);

    /** Returns the zero-based index of the find string in the string list, or -1 if the string
        list does not contain the find string. */
    static int indexOf(const vector<string>& list, string find);

    // ================== Transforms ==================

    /** Replaces all occurences in the text string of 'find' by 'replace', and returns the result.
        */
    static string replace(string text, string find, string replace);

    /** Compresses consecutive white space (space, tab and newline characters) in the text string
        to a single space, removes any white space at the start and at the end, and returns the
        result. */
    static string squeeze(string text);

    /** Pads the text string with the specified pad character at the left until it has the
        specified length (so that the string is right-aligned), and returns the result. The default
        pad character is a space. If the length of the input string is equal to or larger than the
        specified length, the input string is returned unchanged. */
    static string padLeft(string text, size_t length, char pad = ' ');

    /** Pads the text string with the specified pad character at the right until it has the
        specified length (so that the string is left-aligned), and returns the result. The default
        pad character is a space. If the length of the input string is equal to or larger than the
        specified length, the input string is returned unchanged. */
    static string padRight(string text, size_t length, char pad = ' ');

    /** Converts upper-case characters in the text string to lower case, and returns the result. */
    static string toLower(string text);

    /** Converts lower-case characters in the text string to upper case, and returns the result. */
    static string toUpper(string text);

    /** Converts the first character in the string to upper-case, if applicable, and returns the
        result. */
    static string toUpperFirst(string text);

    /** Splits the text string into substrings wherever the separator string occurs, and returns
        the list of those strings. If the separator does not match anywhere, the function returns a
        single-element list containing the complete text string. */
    static vector<string> split(string text, string separator);

    /** Joins the list of strings by inserting the separator between consecutive non-empty strings,
        and returns the resulting single string. */
    static string join(const vector<string>& list, string separator);

    // ================== Filename/path related functions ==================

    /** Joins the two specified path segments into a single path in a meaningful way, and returns
        the result. Specifically, if both segments are nonempty, they are joined such that there is
        exactly one (forward or backward) slash in between. If either (or both) segments are empty,
        no slash is added. */
    static string joinPaths(string segment1, string segment2);

    /** Returns true if the specified path string represents an absolute path; false if it
        represents a relative path. */
    static bool isAbsolutePath(string filepath);

    /** Returns the rightmost path segment in the specified path (after removing a trailing slash or
        backslash). */
    static string filename(string filepath);

    /** Returns the filename proper, without extension, in the specified path (after removing a
        trailing slash or backslash). */
    static string filenameBase(string filepath);

    /** Removes the rightmost path segment from the specified path (leaving a trailing slash or
        backslash in place) and returns the result. If the input string does not contain at least
        two path segments (counting the empty root segment), the empty string is returned. */
    static string dirPath(string filepath);

    /** Adds the specified filename extension to the specified filename (or file path), unless the
        filename already ends with the extension. */
    static string addExtension(string filename, string extension);

    // ================== Conversions ==================

    /** Returns true if the specified string is non-empty and contains a valid representation of a
        boolean. Otherwise returns false. After converting the specified string to lowercase, the
        following contents are considered valid representations: "true", "t", "yes", "y", "1" (for
        boolean true) and "false", "f", "no", "n", "0" (for boolean false). */
    static bool isValidBool(string value);

    /** Returns the boolean value represented by the specified string, or false if the string is
        empty or contains an invalid representation. See isValidBool() for the valid
        representations. */
    static bool toBool(string value);

    /** Returns the string "true" or "false" depending on the specified boolean value. */
    static string toString(bool value);

    /** Returns true if the specified string is non-empty and contains a valid decimal string
        representation of an integer that fits in 32 bits (signed). Otherwise returns false. */
    static bool isValidInt(string value);

    /** Returns the integer value represented by the specified string, or zero if the string is
        empty or contains an invalid representation. See isValidInt() for the valid
        representations. */
    static int toInt(string value);

    /** Returns a string representation of the specified integer value. */
    static string toString(int value);

    /** Returns true if the specified string is non-empty and contains a valid string
        representation of a floating point number that can be represented as a double. Otherwise
        returns false. */
    static bool isValidDouble(string value);

    /** Returns the floating point value represented by the specified string, or zero if the string
        is empty or contains an invalid representation. If the optional argument \em ok is provided
        and not null, the flag is set to true if the conversion succeeded, and to false if not. */
    static double toDouble(string value, bool* ok = nullptr);

    /** Returns a smart string representation of the specified floating point value with up to 10
        significant digits. The function automatically selects fixed point or scientific notation,
        and removes leading and trailing zeroes where applicable, even in the exponent. */
    static string toString(double value);

    /** Returns a string representation of the specified floating point value with the given
        format, precision, and minimum field width.

        The format is 'd' for integer notation (although the value is passed as a double), 'f' for
        fixed point notation, 'e' for scientific notation, and 'g' for the most concise 'f' or 'e'.
        For the 'd' format, the precision is not used and can be omitted (unless the minimum field
        width is specified as well). For the 'f' and 'e' formats, the precision represents the
        number of digits after the decimal point. For the 'g' format, the precision represents the
        maximum number of significant digits (trailing zeroes are omitted). If the format specifier
        is not one of the above, the function silently uses 'd' as a fallback. The precision has a
        default value of 6 and is limited to a maximum of 18 and a minimum of zero.

        If the formatted string is shorter than the specified minimum field width, it is padded
        with the specified pad character so that it becomes right-aligned in the field. The default
        minimum field width is one, which means that no padding occurs. The default pad character
        is a space. The formatted string is never truncated, regardless of the specified width. */
    static string toString(double value, char format, int precision = 6, int width = 1, char pad = ' ');

    /** Returns a user-friendly string representation of the specified memory size value with 3
        significant digits and the appropriate units (KB, MB, GB or TB). */
    static string toMemSizeString(size_t value);
};

////////////////////////////////////////////////////////////////////

#endif
