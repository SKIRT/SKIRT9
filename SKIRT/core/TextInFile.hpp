/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef TEXTINFILE_HPP
#define TEXTINFILE_HPP

#include "Array.hpp"
#include "CompileTimeUtils.hpp"
#include "StoredColumns.hpp"
#include <fstream>
class Log;
class SimulationItem;
class Units;

////////////////////////////////////////////////////////////////////

/** This class allows reading a table of floating point values from a column text file with support
    for unit conversions. The filename is specified in the constructor, and information about the
    columns expected in the file is specified by calling the addColumn() function for each column.

    When reading a column text file, empty lines and lines starting with # are ignored, except for
    the unit specification header lines described below. Other lines must contain a predefined
    number of floating point values separated by whitespace (spaces or tabs). The number of
    expected values is determined for each input file under program control. It is an error for a
    line to contain fewer values than expected (or for any of the values to be improperly
    formatted). On the other hand, any additional information on a line beyond the last expected
    value is ignored.

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
    decribed above), the default units provided by the program are used.

    Stored columns file format
    --------------------------

    Reading large text files can be very slow. Therefore, in certain cases, this class supports the
    possibility to provide the input file in a much faster SKIRT-specific binary file format, i.e.
    the stored columns format described in the StoredColumns class. This file format has the
    following limitations:

    - being a binary format, the files are not human readable.
    - a (binary) file header with column names and unit strings must always be included.
    - the column names and unit strings have a maximum length of 8 characters each.
    - there is no support for non-leaf rows, so the format cannot be used for representing
      adaptive mesh data.

    If the input file provided to the TextInFile constructor has the \c .scol filename extension,
    the implementation automatically switches to reading the binary file format instead of the
    regular column text format. This is fully transparent to the caller of the TextInFile class.
*/
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
        file is successfully opened; (4) \em resource indicates whether the file is located in the
        resources directory (true) or in the regular user input file directory (false, the default
        value).

        If the specified file has the \c .scol filename extension, the implementation automatically
        switches to reading the binary SKIRT column file format instead of the regular column text
        format. For more information, see the class header. */
    TextInFile(const SimulationItem* item, string filename, string description, bool resource = false);

    /** This function closes the file if it was not already closed. It is best to call close() or
        allow the object to go out of scope before logging other messages or starting another
        significant chunk of work. */
    void close();

    /** The destructor calls the close() function. It is best to call close() or allow the object
        to go out of scope before logging other messages or starting another significant chunk of
        work. */
    ~TextInFile();

    //====================== Other functions =======================

    /** This function specifies a mapping (defined by the \em columns argument, as described below)
        between the "physical" columns in the file (defined by the column information in the file
        header) and the "logical" columns requested by the program (defined by repeated calls to
        the addColumn() function). This function can be called with a non-empty \em columns string
        at most once for each file, and such invocation should occur \em before the first
        invocation of the addColumn() function. Calling this function with an empty \em columns
        string is equivalent to not calling it at all.

        If the \em columns string is non-empty, it is interpreted as a comma-separated sequence of
        logical column names. Within each column name, consecutive white space characters are
        replaced by a single space, and white space at the start and at the end is removed. The
        following rules then apply:

        - The input file must contain valid column information in the file header, as described in
        the header of this class.

        - The number of logical column names must match (or exceed) the number of subsequent
        invocations of the addColumn() function.

        - Each logical column name must be equal to exactly one of the file column descriptions,
        unambiguously identifying a particular physical column.

        - Two logical columns cannot identify the same physical column, i.e. a physical column can
        map to at most one logical column.

        - It is allowed for a file to contain physical columns that do not map to a logical column.

        These rules define a mapping between the physical file column ordering and the logical
        column ordering defined by the \em columns string. Once this mapping has been established,
        the program only sees the logical ordering. In other words, the subsequent calls to the
        addColumn() function are matched to the corresponding logical columns, and the readXxx()
        functions retrieve logical columns only. */
    void useColumns(string columns);

    /** This function (virtually) adds a new column to the text file, characterized by the given
        description and unit information. The \em description argument is used only for logging
        purposes. The \em quantity argument specifies the physical quantity represented by the
        column. It must match one of the quantity strings supported by the Units system, or one of
        the special quantity strings recognized by this class (see below). The \em defaultUnit
        argument specifies the default unit string, which is used in case the input file does not
        contain column information.

        In addition to the quantity strings supported by the Units system, this function supports
        the following special quantity strings.

        - The empty string (the default argument value): indicates a dimensionless quantity; the
        default unit must be the empty string as well.

        - The string "specific": indicates a quantity that represents a specific luminosity per
        unit of wavelength, frequency or energy with arbitrary scaling because the values will be
        normalized by the client after being read. The function determines the unit style (per
        wavelength, frequency or energy) based on the units given in the file header or the default
        units. The values are always converted to "per wavelength" style assuming a wavelength
        given by the value of the first preceding column described as "wavelength". However, the
        values will remain scaled with some arbitary wavelength-independent constant.

        The function looks for and, if present, reads the header information line corresponding to
        this column. The unit information from the header is stored with the information provided
        by the function arguments for later use. */
    void addColumn(string description, string quantity = string(), string defaultUnit = string());

    /** This function reads the next row from a column text file and stores the resulting values in
        the array passed to the function by reference. The function first skips empty lines and
        lines starting with a hash character, and then reads a single text line containing data
        values separated by white space.

        The number of expected values corresponds to the number of columns in the file, which is
        determined by repeated calls to the addColumn() function. If the data line contains fewer
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

    /** This is a specialty function intended for use by the AdaptiveMeshSnapshot class when
        importing an adaptive mesh text column file. The function attempts to read a line
        containing a nonleaf node specification. Such a line starts with an exclamation mark, which
        must be followed by three integers (one subdivision specifier for each dimension).

        If the next line (after skipping comments and empty lines) starts with an exclamation mark,
        the function processes it as a nonleaf node specification. If this is successful, the
        function stores the parsed specifiers in the arguments and returns true. If the line cannot
        be parsed, a fatal error is thrown.

        If the next line (after skipping comments and empty lines) does not start with an
        exclamation mark, the contents of the function arguments is undefined and the function
        returns false. In this case, the function has not consumed any information other than
        comments and white space. The file cursor is left just before the next regular line (i.e. a
        line not starting with an exclamation mark), or at the end of the file. */
    bool readNonLeaf(int& nx, int& ny, int& nz);

    /** This variadic template function reads the next row from a column text file and stores the
        resulting values in the variables passed to the function by reference. For example:

        \verbatim
        double a,b,c,d;
        bool success = in.readRow(a,b,c,d);  // reads a row from a file with 4 columns
        \endverbatim

        This function behaves just like the readRow(Array&) version. The number of arguments
        must match the number of columns in the file. */
    template<typename... Values, typename = std::enable_if_t<CompileTimeUtils::isFloatArgList<Values...>()>>
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
    template<typename... Columns, typename = std::enable_if_t<sizeof...(Columns) != 0>>
    void readAllColumns(Columns&... columns)
    {
        auto result = readAllColumns();
        assignColumns(0, result, columns...);
    }

    //======================== Private helpers for column info handling ========================

private:
    /** This function returns the zero-based index of the column that has a file info description
        equal to the given name, or an error value if there is no such column or if there are
        multiple such columns. */
    size_t indexForName(string name) const;

    /** This function returns the index of the first column that is described as "wavelength", or
        the error value if there is no such column. */
    size_t waveIndexForSpecificQuantity() const;

    //======================== Private helpers for reading ========================

private:
    // recursively assign values from Array to double& arguments; used in variadic readRow()
    template<typename... Values>
    static inline void assignValues(size_t index, const Array& result, double& value, Values&... values)
    {
        value = result[index];
        assignValues(index + 1, result, values...);
    }
    static inline void assignValues(size_t /*index*/, const Array& /*result*/) {}

    // recursively assign columns from vector to Array& arguments; used in variadic readAllColumns()
    template<typename... Columns>
    static inline void assignColumns(size_t index, vector<Array>& result, Array& column, Columns&... columns)
    {
        column = std::move(result[index]);
        assignColumns(index + 1, result, columns...);
    }
    static inline void assignColumns(size_t /*index*/, vector<Array>& /*result*/) {}

    //======================== Data Members ========================

private:
    Units* _units{nullptr};  // the units system
    Log* _log{nullptr};      // the logger
    std::ifstream _in;       // the text input stream, if any
    StoredColumns _scol;     // the binary input file, if any

    // private type to store column info
    class ColumnInfo
    {
    public:
        size_t physColIndex{0};  // one-based physical index of this column in the file
        string title;            // description specified in the file, used to remap columns
        string description;      // official description provided by the program
        string quantity;         // quantity, provided by the program
        string unit;             // unit, provided by the program or specified in the file
        double convFactor{1.};   // unit conversion factor from input to internal
        double convPower{1.};    // unit conversion power (exponent) from input to internal
        int waveExponent{0};     // wavelength exponent for converting "specific" quantities
        size_t waveIndex{0};     // zero-based logical index of wavelength column for converting "specific" quantities
    };

    bool _hasTextOpen{false};    // true if a regular column text format file is currently open
    bool _hasBinaryOpen{false};  // true if a binary stored column format file is currently open
    bool _isResource{false};     // true if the file is a resource (as opposed to a user input file)
    bool _hasFileInfo{false};    // becomes true if the file has column header info
    bool _hasProgInfo{false};    // becomes true if the program has added at least one column

    vector<ColumnInfo> _colv;  // info for each column, derived from file info and/or program info
    size_t _numLogCols{0};     // number of logical columns, or number of program columns added so far

    vector<size_t> _logColIndices;  // zero-based index into _colv for each physical column to be read
};

////////////////////////////////////////////////////////////////////

#endif
