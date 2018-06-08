/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "TextInFile.hpp"
#include "FatalError.hpp"
#include "FilePaths.hpp"
#include "Log.hpp"
#include "System.hpp"
#include "Units.hpp"
#include <sstream>

////////////////////////////////////////////////////////////////////

TextInFile::TextInFile(const SimulationItem* item, string filename, string description)
{
    // open the file
    string filepath = item->find<FilePaths>()->input(filename);
    _in = System::ifstream(filepath);
    if (!_in) throw FATALERROR("Could not open the " + description + " text file " + filepath);

    // remember the units system
    _units = item->find<Units>();

    // remember the logger and the message to be issued upon closing
    _log = item->find<Log>();
    _message = item->typeAndName() + " read " + description + " from text file " + filepath + "...";
}

////////////////////////////////////////////////////////////////////

void TextInFile::close()
{
    if (_in.is_open())
    {
        _in.close();
        _log->info(_message);
    }
}

////////////////////////////////////////////////////////////////////

TextInFile::~TextInFile()
{
    close();
}

////////////////////////////////////////////////////////////////////

void TextInFile::addColumn(string description, string quantity, string defaultUnit,
                           bool isGhostColumn, double ghostValue)
{
    _colv.emplace_back(description, quantity, defaultUnit, isGhostColumn,  ghostValue);

    if (quantity=="specific") _colv.back().quantity = "wavelengthmonluminosity";

    // TODO: parse column info in input file and update column info data structures
    // TODO: cache unit conversion factor in column info?
    // TODO: support dimensionless and specific

    /*
    string line;
    while(_in.peek() == '#')
    {
        getline(_in,line);
        if (line.find(find) != string::npos) return line;
    }
    return string();
    */
}

////////////////////////////////////////////////////////////////////

bool TextInFile::readRow(Array& values)
{
    size_t ncols = _colv.size();
    if (!ncols) throw FATALERROR("No columns were declared for column text file");

    // read new line until it is non-empty and non-comment
    string line;
    while (_in.good())
    {
        getline(_in,line);
        auto pos = line.find_first_not_of(" \t");
        if (pos!=string::npos && line[pos]!='#')
        {
            // resize result array if needed (we don't need it to be cleared)
            if (values.size() != ncols) values.resize(ncols);

            // convert values from line and store them in result array
            std::stringstream linestream(line);
            for (size_t i=0; i<ncols; ++i)
            {
                // substitute ghost value if requested for this column
                if (_colv[i].isGhost) values[i] = _colv[i].ghostValue;
                else
                {
                    // convert the value to floating point
                    if (linestream.eof()) throw FATALERROR("One or more required value(s) on text line are missing");
                    double value;
                    linestream >> value;
                    if (linestream.fail()) throw FATALERROR("Input text is not formatted as a floating point number");

                    // convert from input units to internal units
                    // TODO: support dimensionless and specific
                    values[i] = _units->in(_colv[i].quantity, _colv[i].unit, value);
                }
            }
            return true;
        }
    }

    // end of file was reached
    return false;
}

////////////////////////////////////////////////////////////////////

vector<Array> TextInFile::readAllRows()
{
    vector<Array> rows;
    while (true)
    {
        rows.emplace_back();        // add a default-constructed array to the vector
        if (!readRow(rows.back()))  // read next line's values into that array
        {
            rows.pop_back();        // at the end, remove the extraneous array
            break;
        }
    }
    return rows;
}

////////////////////////////////////////////////////////////////////

vector<Array> TextInFile::readAllColumns()
{
    // read the remainder of the file into rows
    const vector<Array>& rows = readAllRows();
    size_t nrows = rows.size();
    size_t ncols = _colv.size();

    // transpose the result into columns
    vector<Array> columns(ncols, Array(nrows));
    for (size_t c=0; c!=ncols; ++c)
        for (size_t r=0; r!=nrows; ++r)
            columns[c][r] = rows[r][c];

    return columns;
}


////////////////////////////////////////////////////////////////////
