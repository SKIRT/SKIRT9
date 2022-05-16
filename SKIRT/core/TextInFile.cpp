/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "TextInFile.hpp"
#include "FatalError.hpp"
#include "FilePaths.hpp"
#include "Log.hpp"
#include "StringUtils.hpp"
#include "System.hpp"
#include "Units.hpp"
#include <exception>
#include <regex>
#include <sstream>

////////////////////////////////////////////////////////////////////

namespace
{
    // This function looks for the next header line that conforms to the required structured syntax. If such a line
    // is found, the column index, description and unit string are stored in the arguments and true is returned.
    // If no such header line is found, the function consumes the complete header and returns false.
    bool getNextInfoLine(std::ifstream& in, size_t& colIndex, string& description, string& unit)
    {
        // continue reading until a conforming header line is found or until the complete header has been consumed
        while (true)
        {
            // consume whitespace characters but nothing else
            while (true)
            {
                auto ch = in.peek();
                if (ch != ' ' && ch != '\t' && ch != '\n' && ch != '\r') break;
                in.get();
            }

            // if the first non-whitespace character is not a hash character, there is no header line
            if (in.peek() != '#') return false;

            // read the header line
            string line;
            getline(in, line);

            // if the line conforms to the required syntax, return the extracted information
            static const std::regex syntax("#\\s*column\\s*(\\d+)\\s*:\\s*([^()]*)\\(\\s*([a-zA-Z0-9/]*)\\s*\\)\\s*",
                                           std::regex::icase);
            std::smatch matches;
            if (std::regex_match(line, matches, syntax) && matches.size() == 4)
            {
                colIndex = std::stoul(matches[1].str());
                description = StringUtils::squeeze(matches[2].str());
                unit = matches[3].str();
                return true;
            }
        }
    }
}

////////////////////////////////////////////////////////////////////

TextInFile::TextInFile(const SimulationItem* item, string filename, string description, bool resource)
{
    // remember the units system and the logger
    _units = item->find<Units>();
    _log = item->find<Log>();

    // get the full path for the resource or for the input file
    _isResource = resource;
    string filepath = resource ? FilePaths::resource(filename) : item->find<FilePaths>()->input(filename);

    // if the file is in SKIRT stored column format, open a binary file
    if (StringUtils::endsWith(filepath, ".scol"))
    {
        _scol.open(filepath);

        // log "reading file" message
        _log->info(item->typeAndName() + " reads " + description + " from binary file " + filepath + "...");

        // read the header information into a list of ColumnInfo records
        auto columnNames = _scol.columnNames();
        auto columnUnits = _scol.columnUnits();
        for (size_t index = 0; index != columnNames.size(); ++index)
        {
            // add a default-constructed ColumnInfo record to the list
            _colv.emplace_back();

            // remember the name and the units specified in the file
            _colv.back().physColIndex = index + 1;  // one-based column index
            _colv.back().unit = columnUnits[index];
            _colv.back().title = columnNames[index];
        }
        _hasFileInfo = true;
        _hasBinaryOpen = true;
    }

    // otherwise open a regular text column file
    else
    {
        _in = System::ifstream(filepath);
        if (!_in) throw FATALERROR("Could not open the " + description + " text file " + filepath);

        // log "reading file" message
        _log->info(item->typeAndName() + " reads " + description + " from text file " + filepath + "...");

        // read any structured header lines into a list of ColumnInfo records
        size_t index;  // one-based column index obtained from file info
        string title;
        string unit;
        while (getNextInfoLine(_in, index, title, unit))
        {
            // add a default-constructed ColumnInfo record to the list
            _colv.emplace_back();
            if (index != _colv.size())
                throw FATALERROR("Incorrect column index in file header for column " + std::to_string(_colv.size()));

            // remember the description and the units specified in the file
            _colv.back().physColIndex = index;
            _colv.back().unit = unit;
            _colv.back().title = title;
        }
        _hasFileInfo = !_colv.empty();
        _hasTextOpen = true;
    }
}

////////////////////////////////////////////////////////////////////

void TextInFile::close()
{
    if (_hasTextOpen) _in.close();
    if (_hasBinaryOpen) _scol.close();

    if (_hasTextOpen || _hasBinaryOpen)
    {
        _hasTextOpen = false;
        _hasBinaryOpen = false;

        // log "done" message, except if this is a resource file or after an exception has been thrown
        if (!_isResource && !std::uncaught_exception()) _log->info("Done reading");
    }
}

////////////////////////////////////////////////////////////////////

TextInFile::~TextInFile()
{
    close();
}

////////////////////////////////////////////////////////////////////

namespace
{
    // Error return values for the functions below
    const int ERROR_NO_EXPON = 999;
    const size_t ERROR_NO_INDEX = 999999;
    const size_t ERROR_AM_INDEX = 999998;

    // This function returns the wavelength exponent needed to convert a per wavelength/frequency/energy
    // quantity to internal (per wavelength) style, given the input units, or the error value if the
    // given units are not supported by any of the relevant quantities.
    int waveExponentForSpecificQuantity(Units* unitSystem, string unitString)
    {
        // a list of known per wavelength / per frequency quantities and the corresponding exponents
        static const vector<string> specificQuantities(
            {"wavelengthmonluminosity", "wavelengthfluxdensity", "wavelengthsurfacebrightness", "neutralmonluminosity",
             "neutralfluxdensity", "neutralsurfacebrightness", "frequencymonluminosity", "frequencyfluxdensity",
             "frequencysurfacebrightness", "energymonluminosity", "energyfluxdensity", "energysurfacebrightness"});
        static const vector<int> specificExponents({0, 0, 0, -1, -1, -1, -2, -2, -2, -3, -3, -3});

        // loop over the list
        for (size_t q = 0; q != specificQuantities.size(); ++q)
        {
            // if this quantity supports the given unit, return the corresponding exponent
            if (unitSystem->has(specificQuantities[q], unitString)) return specificExponents[q];
        }
        return ERROR_NO_EXPON;
    }
}

////////////////////////////////////////////////////////////////////

size_t TextInFile::indexForName(std::string name) const
{
    size_t result = ERROR_NO_INDEX;
    size_t index = 0;
    for (const ColumnInfo& col : _colv)
    {
        if (col.title == name)
        {
            if (result != ERROR_NO_INDEX) return ERROR_AM_INDEX;
            result = index;
        }
        index++;
    }
    return result;
}

////////////////////////////////////////////////////////////////////

size_t TextInFile::waveIndexForSpecificQuantity() const
{
    size_t index = 0;
    for (const ColumnInfo& col : _colv)
    {
        if (col.description == "wavelength") return index;
        index++;
    }
    return ERROR_NO_INDEX;
}

////////////////////////////////////////////////////////////////////

void TextInFile::useColumns(string columns)
{
    // empty columns string behaves as if we were never called at all
    columns = StringUtils::squeeze(columns);
    if (columns.empty()) return;

    // verify that program columns have not yet been added
    if (_hasProgInfo) throw FATALERROR("Program columns were declared before requesting column remapping");

    // verify that file contains column info
    if (!_hasFileInfo) throw FATALERROR("Requesting logical columns but there is no column info in file header");

    // establish the logical column info list
    vector<ColumnInfo> newcolv;
    for (string name : StringUtils::split(columns, ","))
    {
        string sname = StringUtils::squeeze(name);
        size_t index = indexForName(sname);
        if (index == ERROR_NO_INDEX)
            throw FATALERROR("No column description in file header matches logical name '" + sname + "'");
        if (index == ERROR_AM_INDEX)
            throw FATALERROR("Multiple column descriptions in file header match logical name '" + sname + "'");

        newcolv.emplace_back(_colv[index]);
    }

    // replace the column info list
    _colv = newcolv;
}

////////////////////////////////////////////////////////////////////

void TextInFile::addColumn(string description, string quantity, string defaultUnit)
{
    // if the file has no header info at all, add a default record for this column
    if (!_hasFileInfo)
    {
        _colv.emplace_back();
        _colv.back().physColIndex = _numLogCols + 1;
        _colv.back().unit = defaultUnit;
    }
    // otherwise verify that there is a column specification to match this program column
    else
    {
        if (_numLogCols + 1 > _colv.size())
            throw FATALERROR("No column info in file header for column " + std::to_string(_numLogCols + 1));
    }
    _hasProgInfo = true;

    // get a writable reference to the column record being handled, and increment the program column index
    ColumnInfo& col = _colv[_numLogCols++];

    // store the programmatically provided information in the record (unit is already stored)
    col.description = description;
    col.quantity = quantity;

    // verify units and determine conversion factor for this column
    if (col.quantity.empty())  // dimensionless quantity
    {
        if (!col.unit.empty() && col.unit != "1")
            throw FATALERROR("Invalid units (" + col.unit + ") for dimensionless quantity in column "
                             + std::to_string(_numLogCols));
        col.unit = "1";
    }
    else if (col.quantity == "specific")  // arbitrarily scaled value per wavelength or per frequency
    {
        col.waveExponent = waveExponentForSpecificQuantity(_units, col.unit);
        if (col.waveExponent == ERROR_NO_EXPON)
            throw FATALERROR("Invalid units (" + col.unit + ") for specific quantity '" + col.quantity + "' in column "
                             + std::to_string(_numLogCols));
        if (col.waveExponent)
        {
            col.waveIndex = waveIndexForSpecificQuantity();
            if (col.waveIndex == ERROR_NO_INDEX)
                throw FATALERROR("No preceding wavelength column for specific quantity '" + col.quantity
                                 + "' in column " + std::to_string(_numLogCols));
        }
    }
    else
    {
        if (!_units->has(col.quantity, col.unit))
            throw FATALERROR("Invalid units (" + col.unit + ") for quantity '" + col.quantity + "' in column "
                             + std::to_string(_numLogCols));
        double offset;  // all SKIRT units have a zero offset
        std::tie(col.convFactor, col.convPower, offset) = _units->def(col.quantity, col.unit);
    }

    // add the physical to logical column mapping for this column
    if (_logColIndices.size() < col.physColIndex) _logColIndices.resize(col.physColIndex, ERROR_NO_INDEX);
    if (_logColIndices[col.physColIndex - 1] != ERROR_NO_INDEX)
        throw FATALERROR("Multiple logical columns (" + std::to_string(_logColIndices[col.physColIndex - 1] + 1) + ","
                         + std::to_string(_numLogCols) + ") map to the same physical column ("
                         + std::to_string(col.physColIndex) + ")");
    _logColIndices[col.physColIndex - 1] = _numLogCols - 1;

    // for regular user input files, log column information
    if (!_isResource)
    {
        string message = "  Column " + std::to_string(_numLogCols) + ": " + col.description + " (" + col.unit + ")";
        if (!col.title.empty())
        {
            message += " <-- ";
            if (col.physColIndex != _numLogCols) message += "column " + std::to_string(col.physColIndex) + ": ";
            message += col.title;
        }
        _log->info(message);
    }
}

////////////////////////////////////////////////////////////////////

bool TextInFile::readRow(Array& values)
{
    if (!_hasProgInfo) throw FATALERROR("No columns were declared for column text file");

    // read next row in text file
    if (_hasTextOpen)
    {
        // read new line until it is non-empty and non-comment
        string line;
        while (_in.good())
        {
            getline(_in, line);
            auto pos = line.find_first_not_of(" \t");
            if (pos != string::npos && line[pos] != '#')
            {
                // resize result array if needed (we don't need it to be cleared)
                if (values.size() != _numLogCols) values.resize(_numLogCols);

                // convert values from line and store them in result array
                std::stringstream linestream(line);
                for (size_t i : _logColIndices)  // i: zero-based logical index
                {
                    if (linestream.eof()) throw FATALERROR("One or more required value(s) on text line are missing");

                    // read the value as floating point
                    double value;
                    linestream >> value;
                    if (linestream.fail())
                    {
                        // some compilers/libraries do not support reading NaN or Inf values, while others do;
                        // we here provide backup support for NaN values (infinities are more complex because of
                        // the need for handling the negative sign which may be already consumed by the stream)
                        linestream.clear();
                        string offending;
                        linestream >> offending;
                        if (StringUtils::toLower(offending) != "nan")
                            throw FATALERROR("Input text is not formatted as a floating point number: " + offending);
                        value = std::numeric_limits<double>::quiet_NaN();
                    }

                    // if mapped to a logical column, convert from input units to internal units, and store the result
                    if (i != ERROR_NO_INDEX)
                    {
                        const ColumnInfo& col = _colv[i];
                        if (col.convPower != 1.) value = pow(value, col.convPower);
                        value *= (col.waveExponent ? pow(values[col.waveIndex], col.waveExponent) : col.convFactor);
                        values[i] = value;
                    }
                }
                return true;
            }
        }
    }

    // read next row in binary file
    else if (_hasBinaryOpen)
    {
        const double* row = _scol.nextRow();
        if (row)
        {
            // resize result array if needed (we don't need it to be cleared)
            if (values.size() != _numLogCols) values.resize(_numLogCols);

            // copy values from row into result array
            for (size_t i : _logColIndices)  // i: zero-based logical index
            {
                // get the value
                double value = *row++;

                // if mapped to a logical column, convert from input units to internal units, and store the result
                if (i != ERROR_NO_INDEX)
                {
                    const ColumnInfo& col = _colv[i];
                    if (col.convPower != 1.) value = pow(value, col.convPower);
                    value *= (col.waveExponent ? pow(values[col.waveIndex], col.waveExponent) : col.convFactor);
                    values[i] = value;
                }
            }
            return true;
        }
    }

    // end of file was reached or no file is open
    return false;
}

////////////////////////////////////////////////////////////////////

bool TextInFile::readNonLeaf(int& nx, int& ny, int& nz)
{
    // read next non-leaf row in text file
    if (_hasTextOpen)
    {
        string line;

        while (true)
        {
            int c = _in.peek();

            // skip comments line
            if (c == '#')
            {
                getline(_in, line);
            }

            // process nonleaf line
            else if (c == '!')
            {
                _in.get();  // skip exclamation mark
                getline(_in, line);

                // convert nx,ny,nz values from line and store them in output arguments
                std::stringstream linestream(line);
                linestream >> nx >> ny >> nz;
                if (linestream.fail())
                    throw FATALERROR("Nonleaf subdivision specifiers are missing or not formatted as integers");

                return true;
            }

            // eat leading white space and empty lines
            else if (c == ' ' || c == '\t' || c == '\n' || c == '\r')
            {
                _in.get();
            }

            // signal not a nonleaf line
            else
            {
                return false;
            }
        }
    }

    // unsupported binary file format, or no file is open
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
            rows.pop_back();  // at the end, remove the extraneous array
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
    for (size_t c = 0; c != ncols; ++c)
        for (size_t r = 0; r != nrows; ++r) columns[c][r] = rows[r][c];

    return columns;
}

////////////////////////////////////////////////////////////////////
