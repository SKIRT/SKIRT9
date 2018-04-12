/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef TEXTINFILE_HPP
#define TEXTINFILE_HPP

#include "Array.hpp"
#include <fstream>
class SimulationItem;

////////////////////////////////////////////////////////////////////

/** This class allows reading floating point values from the input text file specified in the
    constructor. The values should be organized in columns, forming a table. An informational
    message is logged when the file is opened, and the file is automatically closed when the object
    is destructed. */
class TextInFile
{
    //=============== Construction - Destruction  ==================

public:
    /** The constructor opens the specified file for reading, and logs a message when successful.
        If the file can't be opened, a FatalError is thrown. The constructor takes several
        arguments: (1) \em item specifies a simulation item in the hierarchy of the caller (usually
        the caller itself) used to retrieve the input file path and an appropriate logger; (2) \em
        filename specifies the name of the file, including filename extension but excluding path
        and simulation prefix; (3) \em description specifies a description used in the log message
        issued after the file is successfully opened; */
    TextInFile(const SimulationItem* item, string filename, string description);

    //====================== Other functions =======================

    /** This function looks through the header at the current position in the file, and returns the
        first header line containing the specified string, or the empty string if no such line is
        found. In this context, a header is defined as consecutive set of lines starting with a
        hash. A line not starting with a hash (including an empty line) ends the header and will
        not be consumed by this function. */
    string readHeaderLine(string find);

    /** This function reads the next row from a column text file and stores the resulting values in
        the array passed to the function by reference. The function first skips empty lines and
        lines starting with a hash character, and then reads a single text line containing between
        \em ncols - \em noptcols and up to \em ncols values separated by white space. These values
        are converted to floating point and stored in the specified \em values array. If a row was
        successfully read, the \em values array is resized to \em ncols and overwritten with the
        results, and the function returns true. Missing optional values at the end of the row are
        replaced by zeroes. Any additional information on the line beyond the last value is
        ignored. If the end of the file is reached before a row can be read, the function returns
        false and the contents of the \em values array is undefined. If the line contains
        improperly formatted floating point numbers, or if there are less than \em ncols - \em
        noptcols values, the function throws a FatalError (and the contents of the \em values array
        is undefined). */
    bool readRow(Array& values, size_t ncols, size_t noptcols = 0);

    /** This variadic template function reads the next row from a column text file and stores the
        resulting values in the variables passed to the function by reference. For example:

        \verbatim
        double a,b,c,d;
        bool success = in.readRow(1, a,b,c,d);  // reads a line with 4 values,
                                                // the last of which is optional
        \endverbatim

        This function behaves just like the readRow(Array&) version, except that the number of
        requested values is derived from the number of double& arguments, and that the number of
        optional values (\em noptcols) is specified as the first argument. */
    template <typename... Values>
    bool readRow(size_t noptcols, Values&... values);

    /** This function reads all rows from a column text file (from the current position until the
        end of the file), and returns the resulting values as a vector of arrays. For each row,
        this function behaves just like readRow(Array&). */
    std::vector<Array> readAllRows(size_t ncols, size_t noptcols = 0);

private:
    // recursively assign values from Array to double& arguments; used in variadic readRow()
    template <typename... Values>
    void assign(size_t index, const Array& result, double& value, Values&... values);
    inline void assign(size_t index, const Array& result);

    //======================== Data Members ========================

private:
    std::ifstream _in;  // the input stream
};

////////////////////////////////////////////////////////////////////

template <typename... Values>
bool TextInFile::readRow(size_t noptcols, Values&... values)
{
    Array result;
    bool success = readRow(result, sizeof...(Values), noptcols);
    if (success) assign(0, result, values...);
    return success;
}

template <typename... Values>
void TextInFile::assign(size_t index, const Array& result, double& value, Values&... values)
{
    value = result[index];
    assign(index+1, result, values...);
}
void TextInFile::assign(size_t /*index*/, const Array& /*result*/)
{
}

////////////////////////////////////////////////////////////////////

#endif
