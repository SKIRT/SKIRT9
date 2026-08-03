/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef STOREDCOLUMNS_HPP
#define STOREDCOLUMNS_HPP

#include "Array.hpp"

////////////////////////////////////////////////////////////////////

/** An instance of the StoredColumns class provides access to a set of data columns stored in a
    file using the binary SKIRT stored columns format (which is simular to the SKIRT stored table
    format implemented by StoredTable). Stored columns files are intended as a much faster
    replacement for large regular text column data files, without the benefit of being human
    readable. The format does not support non-leaf rows, so it cannot be used for representing
    adaptive mesh data.

    Stored columns file format
    --------------------------

    A stored columns file includes the names and units of the quantities in each column, in
    addition to the tabulated data values. All values are stored as binary data in the form of
    64-bit floating-point numbers. More specifically, a stored columns file is essentially a
    sequence of 8-byte data items. A data item can have one of three types:

        - string: 1 to 8 printable and non-whitespace 7-bit ASCII characters, padded with spaces
          to fill 8 bytes if needed;
        - unsigned integer: 64-bit integer in little-endian byte order;
        - floating point: 64-bit double (IEEE 754) in little-endian byte order.

    The overall layout is as follows:
        - SKIRT name/version tag
        - Endianness tag
        - 0  (to differentiate from stored table format, which stores the nr of axes here)
        - numRows
        - numColumns
        - columnName (x numColumns)
        - columnUnit (x numColumns)
        - value (x numColumns x numRows)
        - end-of-file tag

    The values are ordered so that the column values for a particular row are next to each other.

    The StoredColumns class
    -----------------------

    The default constructor creates an invalid stored columns instance. The alternate constructor
    and the open() function associate a particular stored columns file with the stored columns
    instance. The close() function and the destructor automatically release the file association
    and any related resources. It is allowed to call open() again after close().

    The columnNames() and columnUnits() functions return informaton about the columns, and the
    readRow() function returns the rows one by one, from the start to the end of the file.
*/
class StoredColumns
{
    // ================== Constructing ==================

public:
    /** The default constructor constructs a closed stored columns instance. The user must call the
        open() function to associate the instance with a particular stored columns file. */
    StoredColumns() {}

    /** This alternate constructor constructs a stored columns instance and immediately associates
        a given stored columns file with it by calling the open() function. Refer to the open()
        function for more information. */
    StoredColumns(string filename) { open(filename); }

    /** The destructor releases the association with a stored columns file established by the
        alternate constructor or the open() function, if there is any. */
    ~StoredColumns() { close(); }

    /** This function associates a given stored columns input file with the stored columns
        instance. If such an association already exists, or if the open operations fails, this
        function throws a fatal error. The \em filename argument specifies the absolute file path
        of the input file, including the mandatory ".scol" filename extension. */
    void open(string filepath);

    /** This function releases the association with a stored columns file established by the
        alternate constructor or the open() function, if there is any. After callig this function,
        it is allowed to call open() again. */
    void close();

    // ================== Accessing data ==================

public:
    /** This function returns the names of the columns in the file. The length of the list
        corresponds to the number of columns. If no file is open, the list is enpty. */
    const vector<string>& columnNames() const { return _columnNames; }

    /** This function returns the unit strings for the columns in the file, in the same order as
        the column names returned by columnNames(). If no file is open, the list is enpty. */
    const vector<string>& columnUnits() const { return _columnUnits; }

    /** This function returns a pointer to a list of the column values in the next row, or the null
        pointer if there are no more rows or if no file is open. The returned pointer becomes
        invalid upon the next invocation of this function, or when the file is closed. The column
        values are listed in the order corresponding to the list returned by columnNames(). The
        function returns rows one by one, in the order they are stored in the file. */
    const double* nextRow()
    {
        if (_nextRow == _endRow) return nullptr;
        auto result = _nextRow;
        _nextRow += _numColumns;
        return result;
    }

    // ================== Data members ==================

private:
    string _filePath;                 // the canonical path to the associated stored columns file
    vector<string> _columnNames;      // the column names
    vector<string> _columnUnits;      // the column unit strings
    size_t _numColumns{0};            // step size from one row to the next
    const double* _nextRow{nullptr};  // ptr to the first value in the next row to be returned
    const double* _endRow{nullptr};   // ptr just beyond the last row in the file
};

////////////////////////////////////////////////////////////////////

#endif
