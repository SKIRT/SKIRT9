/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef TEXTINFILE_HPP
#define TEXTINFILE_HPP

#include "Array.hpp"
#include "CompileTimeUtils.hpp"
#include <fstream>
class Log;
class SimulationItem;
class Units;
namespace TextInFile_Private { class ColumnInfo; }

////////////////////////////////////////////////////////////////////

/** This class allows reading a table of floating point values from a column text file with support
    for unit conversions. The filename is specified in the constructor, and information about the
    columns expected in the file is specified by calling the addColumn() function for each column.
    It is recommended that the header of the input file includes column information in the format
    described below. For example, a two-column input file might have the following header:

        # My personal SED written from Python
        # column 1: wavelength (A)
        # column 2: specific luminosity (W/Hz)
        #
        ... data ...

    There must be a header line for each column with the following contents from left to right:
      - the hash character indicating a header line
      - the word "column" (case insensitive)
      - the one-based column number (optional)
      - a colon
      - a description that does not contain parenthesis (optional)
      - a unit string between parenthesis (may be "(1)" or "()" for dimensionless quantities)

    Extra whitespace is allowed between the various fields. The unit string must exactly match one
    of the unit strings supported by SKIRT for the expected quantity. The column info lines must
    occur in the same order as the data columns. Extra header lines are allowed between, before,
    and after the column info lines.

    If there is no column information in the file (i.e. none of the header lines match the syntax
    decribed above), the default units provided by the program are used. */
class TextInFile
{
    //=============== Construction - Destruction  ==================

public:
    /** The constructor opens the specified file for reading; if the file can't be opened, a
        FatalError is thrown. The constructor takes several arguments: (1) \em item specifies a
        simulation item in the hierarchy of the caller (usually the caller itself) used to retrieve
        the input file path and an appropriate logger; (2) \em filename specifies the name of the
        file, including filename extension but excluding path and simulation prefix; (3) \em
        description describes the contents of the file for use in the log message issued after the
        file is successfully opened. */
    TextInFile(const SimulationItem* item, string filename, string description);

    /** This function closes the file if it was not already closed. It is best to call close() or
        allow the object to go out of scope before logging other messages or starting another
        significant chunk of work. */
    void close();

    /** The destructor calls the close() function. It is best to call close() or allow the object
        to go out of scope before logging other messages or starting another significant chunk of
        work. */
    ~TextInFile();

    //====================== Other functions =======================

    /** This function (virtually) adds a new column to the text file, characterized by the given
        description and unit information. The \em description argument is used only for
        logging purposes. The \em quantity argument specifies the physical quantity represented by
        the column. It must match one of the quantity strings supported by the Units system, or one
        of the special quantity strings recognized by this class (see below). The \em defaultUnit
        argument specifies the default unit string, which is used in case the input file does not
        contain column information.

        In addition to the quantity strings supported by the Units system, this function supports
        the following special quantity strings.
           - The empty string (the default argument value): indicates a dimensionless quantity;
             the default unit must be the empty string as well.
           - The string "specific": indicates a quantity that represents a specific luminosity per
             unit of frequency or per unit of wavelength, in arbitrary units (because the values
             will be normalized after being read). The function determines the frequency/wavelength
             flavor based on the units given in the file header or the default units. The values
             are converted to "per wavelength" flavor if needed using the value of the first
             preceding column described as "wavelength". However, the values will remain scaled
             with some arbitary wavelength-independent constant.

        For a regular column, the \em isGhostColumn flag must be omitted or set to false. In this
        case, this function looks for and, if present, reads the header information line
        corresponding to this column. The unit information from the header is stored with the
        information provided by the function arguments for later use.

        On the other hand, if the \em isGhostColumn flag is set to true, the function constructs a
        \em ghost \em column that has no counterpart in the input file. Rather than reading a value
        from the file, the value specified as the \em ghostValue argument is copied into the output
        array instead. This capability is useful to keep the number of values obtained by readRow()
        constant even if some of the columns are sometimes missing from the input file.
    */
    void addColumn(string description, string quantity = string(), string defaultUnit = string(),
                   bool isGhostColumn = false, double ghostValue = 0.);

    /** This function reads the next row from a column text file and stores the resulting values in
        the array passed to the function by reference. The function first skips empty lines and
        lines starting with a hash character, and then reads a single text line containing data
        values separated by white space.

        The number of expected values corresponds to the number of columns in the file, which is
        determined by repeated calls to the addColumn() function. If the data line contains less
        values than expected, or if any of the values is improperly formatted for a floating point
        number, the function throws a FatalError (the size and contents of the \em values array are
        undefined). Any additional information on the line beyond the last expected value is
        ignored.

        If a row was successfully read, the input values are converted from the input units (as
        specified in the file header or using the default given in the addColumn() function) to
        SKIRT-internal units. Finally, the \em values array is set to the appropriate length, the
        converted input values are stored into it in column order, and the function returns true.

        If the end of the file is reached before a row can be read, the function returns false and
        the size and contents of the \em values array are undefined. */
    bool readRow(Array& values);

    /** This variadic template function reads the next row from a column text file and stores the
        resulting values in the variables passed to the function by reference. For example:

        \verbatim
        double a,b,c,d;
        bool success = in.readRow(a,b,c,d);  // reads a row from a file with 4 columns
        \endverbatim

        This function behaves just like the readRow(Array&) version. The number of arguments
        must match the number of columns in the file. */
    template <typename... Values, typename = std::enable_if_t<CompileTimeUtils::isFloatArgList<Values...>()>>
    bool readRow(Values&... values)
    {
        Array result;
        bool success = readRow(result);
        if (success) assignValues(0, result, values...);
        return success;
    }

    /** This function reads all rows from a column text file (from the current position until the
        end of the file), and returns the resulting values as a vector of row arrays. For each row,
        this function behaves just like readRow(Array&). */
    vector<Array> readAllRows();

    /** This function reads all rows from a column text file (from the current position until the
        end of the file), transposes the data repesentation from rows into columns, and returns the
        resulting values as a vector of column arrays. For each row, this function behaves just
        like readRow(Array&). */
    vector<Array> readAllColumns();

    /** This function reads all rows from a column text file (from the current position until the
        end of the file), transposes the data repesentation from rows into columns, and stores the
        resulting column arrays in the variables passed to the function by reference. For each row,
        this function behaves just like readRow(Array&). */
    template <typename... Columns, typename = std::enable_if_t<sizeof...(Columns)!=0>>
    void readAllColumns(Columns&... columns)
    {
        auto result = readAllColumns();
        assignColumns(0, result, columns...);
    }

private:
    // recursively assign values from Array to double& arguments; used in variadic readRow()
    template <typename... Values>
    void assignValues(size_t index, const Array& result, double& value, Values&... values)
    {
        value = result[index];
        assignValues(index+1, result, values...);
    }
    inline void assignValues(size_t /*index*/, const Array& /*result*/) { }

    // recursively assign columns from vector to Array& arguments; used in variadic readAllColumns()
    template <typename... Columns>
    void assignColumns(size_t index, vector<Array>& result, Array& column, Columns&... columns)
    {
        column = std::move(result[index]);
        assignColumns(index+1, result, columns...);
    }
    inline void assignColumns(size_t /*index*/, vector<Array>& /*result*/) { }

    //======================== Data Members ========================

private:
    // data members initialized in the constructor
    std::ifstream _in;      // the input stream
    Units* _units{nullptr}; // the units system
    Log* _log{nullptr};     // the logger

    // data members initialized by repeated calls to addColumn()
    vector<TextInFile_Private::ColumnInfo*> _colv;  // info for each column
    bool _hasHeaderInfo{true};   // true if the input file has structured column info in header
};

////////////////////////////////////////////////////////////////////

#endif
