/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "TextInFile.hpp"
#include "FatalError.hpp"
#include "FilePaths.hpp"
#include "Log.hpp"
#include "System.hpp"
#include <sstream>

////////////////////////////////////////////////////////////////////

TextInFile::TextInFile(const SimulationItem* item, string filename, string description)
{
    // open the file and log a message
    string filepath = item->find<FilePaths>()->input(filename);
    _in = System::ifstream(filepath);
    if (!_in) throw FATALERROR("Could not open the " + description + " data file " + filepath);
    item->find<Log>()->info("Reading " + description + " from file " + filepath + "...");
}

////////////////////////////////////////////////////////////////////

string TextInFile::readHeaderLine(string find)
{
    string line;
    while(_in.peek() == '#')
    {
        getline(_in,line);
        if (line.find(find) != string::npos) return line;
    }
    return string();
}

////////////////////////////////////////////////////////////////////

bool TextInFile::readRow(Array& values, size_t ncols, size_t noptcols)
{
    // read new line until it is non-empty and non-comment
    string line;
    while (_in.good())
    {
        getline(_in,line);
        auto pos = line.find_first_not_of(" \t");
        if (pos!=string::npos && line[pos]!='#')
        {
            // resize and clear result array
            values.resize(ncols);

            // convert values from line and store them in result array
            std::stringstream linestream(line);
            for (size_t i=0; i<ncols; ++i)
            {
                if (linestream.eof())
                {
                    if (i<ncols-noptcols) throw FATALERROR("One or more required value(s) on text line are missing");
                    break;
                }
                double value;
                linestream >> value;
                if (linestream.fail()) throw FATALERROR("Input text is not formatted as a floating point number");
                values[i] = value;
            }
            return true;
        }
    }

    // end of file was reached
    return false;
}

////////////////////////////////////////////////////////////////////

vector<Array> TextInFile::readAllRows(size_t ncols, size_t noptcols)
{
    vector<Array> rows;
    while (true)
    {
        rows.emplace_back();                        // add a default-constructed array to the vector
        if (!readRow(rows.back(), ncols, noptcols)) // read next line's values into that array
        {
            rows.pop_back();                        // at the end, remove the extraneous array
            break;
        }
    }
    return rows;
}

////////////////////////////////////////////////////////////////////
