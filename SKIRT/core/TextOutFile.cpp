/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "TextOutFile.hpp"
#include "FatalError.hpp"
#include "FilePaths.hpp"
#include "Log.hpp"
#include "ProcessManager.hpp"
#include "StringUtils.hpp"
#include "System.hpp"
#include "Units.hpp"

////////////////////////////////////////////////////////////////////

TextOutFile::TextOutFile(const SimulationItem* item, string filename, string description, bool overwrite)
{
    _log = item->find<Log>();
    _filepath = item->find<FilePaths>()->output(filename + ".dat");
    _units = item->find<Units>();

    // Only open the output file if this is the root process
    if (ProcessManager::isRoot())
    {
        _log->info("Writing " + description + " to " + _filepath + "...");
        _out = System::ofstream(_filepath, !overwrite);
    }
}

////////////////////////////////////////////////////////////////////

TextOutFile::~TextOutFile()
{
    // Close the output file, if it was opened by this process
    if (_out.is_open())
    {
        _out.close();
        _log->info("File " + _filepath + " created.");
    }
}

////////////////////////////////////////////////////////////////////

void TextOutFile::addColumn(string description, char format, int precision)
{
    _formats.push_back(format);
    _precisions.push_back(precision);

    writeLine("# column " + std::to_string(++_ncolumns) + ": " + description);
}

////////////////////////////////////////////////////////////////////

void TextOutFile::writeLine(string line)
{
    if (_out.is_open())
    {
        _out << line << std::endl;
    }
}

////////////////////////////////////////////////////////////////////

void TextOutFile::writeRow(vector<double> values)
{
    if (values.size() != _ncolumns) throw FATALERROR("Number of values in row does not match the number of columns");

    string line;
    for (size_t i=0; i<_ncolumns; i++)
    {
        line += (i ? " " : "") + StringUtils::toString(values[i], _formats[i], _precisions[i]);
    }
    writeLine(line);
}

////////////////////////////////////////////////////////////////////
