/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef TEXTOUTFILE_HPP
#define TEXTOUTFILE_HPP

#include "Array.hpp"
#include "CompileTimeUtils.hpp"
#include <array>
#include <fstream>
class Log;
class SimulationItem;
class Units;

////////////////////////////////////////////////////////////////////

/** This class allows writing text to a file specified in the constructor, with explicit support
    for formatting columns of floating point or integer numbers. Text is written per line, by
    calling the writeLine() or writeRow() functions. In a multiprocessing environment, only the
    root process will be allowed to write to the specified file; calls to writeLine() or writeRow()
    performed by other processes will have no effect. */
class TextOutFile
{
    //=============== Construction - Destruction  ==================

public:
    /** The constructor of the TextOutFile class. If the constructor is invoked from the root
        process, the output file is opened and a stream for the file is initialized. For other
        processes, the output stream remains uninitialized and is never used, i.e. calling the
        writeLine() function has no effect. The constructor takes several arguments: (1) \em item
        specifies a simulation item in the hierarchy of the caller (usually the caller itself) used
        to retrieve the output file path and an appropriate logger, and to determine whether this
        is the root process; (2) \em filename specifies the name of the file, excluding path,
        simulation prefix and filename extension; (3) \em description describes the contents of the
        file for use in the log message issued after the file is successfully closed. */
    TextOutFile(const SimulationItem* item, string filename, string description);

    /** In the root process, this function closes the file and logs an informational message, if
        the file was not already closed. It is important to call close() or allow the object to go
        out of scope before logging other messages or starting another significant chunk of work.
        */
    void close();

    /** The destructor calls the close() function. It is important to call close() or allow the
        object to go out of scope before logging other messages or starting another significant
        chunk of work. */
    ~TextOutFile();

    //====================== Other functions =======================

public:
    /** This function writes the specified string to the file as a new line. If the calling process
        is not the root, this function will have no effect. */
    void writeLine(string line);

    /** This function (virtually) adds a new column to the text file, characterized by a certain
        description and formatting. The format is 'd' for integer values, 'e' for scientific notation, 'f'
        for floating point notation and 'g' for the most concise 'e' or 'f'. For the
        'e' and 'f' formats, the precision represents the number of digits after the decimal
        point. For the 'g' format, the precision represents the maximum number of significant
        digits (trailing zeroes are omitted). The description of each column is added to the header of
        the text file, along with the column number. */
    void addColumn(string quantityDescription, string unitDescription = string(), char format = 'e', int precision = 9);

    /** This function writes the specified list of (double) values to the text file, on a single row
        where adjacent values are seperated by a space. The values are formatted according to the
        'format' and 'precision' specified by the addColumn function. If the number of values in the
        list does not match the number of columns, a FatalError is thrown. */
    void writeRow(vector<double> values);

    /** This template function writes the specified list of values to the text file, on a single
        row where adjacent values are seperated by a space. The values are formatted according to
        the 'format' and 'precision' specified by the addColumn function. If the number of values
        in the list does not match the number of columns, a FatalError is thrown. */
    template<typename... Values, typename = std::enable_if_t<CompileTimeUtils::isNumericArgList<Values...>()>>
    void writeRow(Values... values)
    {
        std::array<double, sizeof...(values)> list = {{static_cast<double>(values)...}};
        writeRowPrivate(sizeof...(values), &list[0]);
    }

private:
    /** This function writes the specified list of (double) values to the text file with the same
        semantics as the other writeRow() functions. It is intended for private use from the
        template writeRow() functions. */
    void writeRowPrivate(size_t n, const double* values);

    //======================== Data Members ========================

protected:
    // can be used by subclasses
    Units* _units{nullptr};  // for conversion to output units
    std::ofstream _out;      // the output stream

private:
    // used for column formatting
    size_t _ncolumns{0};
    vector<char> _formats;
    vector<int> _precisions;

    // used when closing
    Log* _log{nullptr};  // the logger
    string _message;     // the message
};

////////////////////////////////////////////////////////////////////

#endif
