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
#include <exception>

////////////////////////////////////////////////////////////////////

TextOutFile::TextOutFile(const SimulationItem* item, string filename, string description)
{
    // Only open the output file if this is the root process
    if (ProcessManager::isRoot())
    {
        // open the file
        string filepath = item->find<FilePaths>()->output(filename + ".dat");
        _out = System::ofstream(filepath);
        if (!_out) throw FATALERROR("Could not open the " + description + " output file " + filepath);

        // remember some pointers
        _log = item->find<Log>();
        _units = item->find<Units>();

        // remember the message to be issued upon closing
        _message = item->typeAndName() + " wrote " + description + " to " + filepath;
    }
}

////////////////////////////////////////////////////////////////////

void TextOutFile::close()
{
    if (_out.is_open())
    {
        _out.close();

        // log success message, except if an exception has been thrown
        if (!std::uncaught_exception()) _log->info(_message);
        ;
    }
}

////////////////////////////////////////////////////////////////////

TextOutFile::~TextOutFile()
{
    close();
}

////////////////////////////////////////////////////////////////////

void TextOutFile::addColumn(string quantityDescription, string unitDescription, char format, int precision)
{
    _formats.push_back(format);
    _precisions.push_back(precision);

    if (unitDescription.empty()) unitDescription = "1";
    writeLine("# column " + std::to_string(++_ncolumns) + ": " + quantityDescription + " (" + unitDescription + ")");
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
    writeRowPrivate(values.size(), &values[0]);
}

////////////////////////////////////////////////////////////////////

void TextOutFile::writeRowPrivate(size_t n, const double* values)
{
    if (n != _ncolumns) throw FATALERROR("Number of values in row does not match the number of columns");

    string line;
    for (size_t i = 0; i < _ncolumns; i++)
    {
        line += (i ? " " : "") + StringUtils::toString(values[i], _formats[i], _precisions[i]);
    }
    writeLine(line);
}

////////////////////////////////////////////////////////////////////
