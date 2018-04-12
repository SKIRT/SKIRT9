/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef TEXTOUTFILE_HPP
#define TEXTOUTFILE_HPP

#include "Array.hpp"
#include "SimulationItem.hpp"
#include <fstream>
class Log;
class Units;

////////////////////////////////////////////////////////////////////

/** This class represents a writable text file that can be initizialized by providing a filename to
    its constructor. Text is written per line, by calling the writeLine() function. The file is
    automatically closed when the object is destructed. In a multiprocessing environment, only the
    root process will be allowed to write to the specified file; calls to writeLine() performed by
    other processes will have no effect. */
class TextOutFile
{
    //=============== Construction - Destruction  ==================

public:
    /** The constructor of the TextOutFile class. If the constructor is invoked from the root
        process, the output stream for the file is initialized and a log message is issued. For
        other processes, this output stream remains uninitialized and is never used, i.e. calling
        the writeLine() function has no effect. The constructor takes several arguments: (1) \em
        item specifies a simulation item in the hierarchy of the caller (usually the caller itself)
        used to retrieve the output file path and an appropriate logger; (2) \em filename specifies
        the name of the file, excluding path, simulation prefix and filename extension; (3) \em
        description specifies a description used in the log message issued after the file was
        successfully written; (4) \em overwrite is an optional flag that can be set to false to
        append the new lines to the existing file (it should be left to the default of true in
        almost all cases). */
    TextOutFile(const SimulationItem* item, string filename, string description, bool overwrite = true);

    /** The destructor of the TextOutFile class. On the root process, the file is closed and a log
        message is issued. */
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
    void addColumn(string description, char format = 'e', int precision = 6);

    /** This function writes the specified list of (double) values to the text file, on a single row
        where adjecent values are seperated by a space. The values are formatted according to the
        'format' and 'precision' specified by the addColumn function. If the number of values in the
        list does not match the number of columns, a FatalError is thrown. */
    void writeRow(vector<double> values);

    //======================== Data Members ========================

protected:
    Log* _log;
    Units* _units;
    std::ofstream _out;  // the output stream

private:
    string _filepath;
    size_t _ncolumns{0};
    vector<char> _formats;
    vector<int> _precisions;
};

////////////////////////////////////////////////////////////////////

#endif
