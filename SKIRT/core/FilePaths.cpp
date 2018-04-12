/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "FilePaths.hpp"
#include "FatalError.hpp"
#include "StringUtils.hpp"
#include "System.hpp"
#include <mutex>

////////////////////////////////////////////////////////////////////

namespace
{
    // flag becomes true if the static paths have been initialized
    std::once_flag _initialized;

    // the static application and resource paths
    string _applicationPath;
    string _resourcePath;

    // relative paths to check for presence of data folder (built-in resources)
    const char* _datpaths[] = { "data", "../../../git/SKIRT/data", "../../../../git/SKIRT/data" };
    const int _Ndatpaths = sizeof(_datpaths) / sizeof(const char*);

    // relative paths to check for presence of extdat folder (external resources)
    const char* _extdatpaths[] = { "extdat", "../../../extdat", "../../../../extdat" };
    const int _Nextdatpaths = sizeof(_extdatpaths) / sizeof(const char*);

    // sets the static application and resource paths, or throws an error if there is a problem
    void setStaticPaths()
    {
        // get the executable path (or the empty string in case of failure)
        string execPath = System::executablePath();
        if (execPath.empty()) throw FATALERROR("Could not determine path to executable");

        // store the location of the executable (i.e. the path to the containing directory)
        _applicationPath = StringUtils::dirPath(execPath);

        // iterate over the relative paths
        for (int i=0; i<_Ndatpaths; i++)
        {
            string test = StringUtils::joinPaths(_applicationPath, _datpaths[i]);
            if (System::isDir(test))
            {
                _resourcePath = System::canonicalPath(test) + "/";
                return;
            }
        }

        // if we reach here, the resource folder wasn't found
        throw FATALERROR("Could not locate 'data' directory relative to '" + _applicationPath + "'");
    }
}

////////////////////////////////////////////////////////////////////

FilePaths::FilePaths(SimulationItem* parent)
{
    parent->addChild(this);
}

////////////////////////////////////////////////////////////////////

void FilePaths::setupSelfBefore()
{
    SimulationItem::setupSelfBefore();
    std::call_once(_initialized, setStaticPaths);
}

////////////////////////////////////////////////////////////////////

void FilePaths::setInputPath(string value)
{
    if (!System::isDir(value)) throw FATALERROR("Input path does not exist or is not a directory: " + value);
    _inputPath = System::canonicalPath(value) + "/";
}

////////////////////////////////////////////////////////////////////

string FilePaths::inputPath() const
{
    return _inputPath;
}

////////////////////////////////////////////////////////////////////

void FilePaths::setOutputPath(string value)
{
    if (!System::isDir(value)) throw FATALERROR("Output path does not exist or is not a directory: " + value);
    _outputPath = System::canonicalPath(value) + "/";
}

////////////////////////////////////////////////////////////////////

string FilePaths::outputPath() const
{
    return _outputPath;
}

////////////////////////////////////////////////////////////////////

void FilePaths::setOutputPrefix(string value)
{
    _outputPrefix = value;
}

////////////////////////////////////////////////////////////////////

string FilePaths::outputPrefix() const
{
    return _outputPrefix;
}

////////////////////////////////////////////////////////////////////

string FilePaths::input(string name) const
{
    return _inputPath + name;
}

////////////////////////////////////////////////////////////////////

string FilePaths::output(string name) const
{
    return _outputPath + _outputPrefix + "_" + name;
}

////////////////////////////////////////////////////////////////////

string FilePaths::resource(string name)
{
    // initialize the static paths if needed
    std::call_once(_initialized, setStaticPaths);
    return _resourcePath + name;
}

////////////////////////////////////////////////////////////////////

namespace
{
    // recursively searches the given directory and its subdirectories for the given file
    // returns the full path if found or the empty string if not
    string findFile(string directory, string filename)
    {
        // search the current level
        if (StringUtils::contains(System::filesInDirectory(directory), filename))
        {
            return StringUtils::joinPaths(directory, filename);
        }

        // search the subdirectories
        for (const string& subdir : System::dirsInDirectory(directory))
        {
            string result = findFile(StringUtils::joinPaths(directory, subdir), filename);
            if (!result.empty()) return result;
        }

        // report failure
        return string();
    }
}

////////////////////////////////////////////////////////////////////

string FilePaths::externalResource(string name)
{
    // initialize the static paths if needed
    std::call_once(_initialized, setStaticPaths);

    // iterate over the relative paths
    for (int i=0; i<_Nextdatpaths; i++)
    {
        // recursively iterate over subdirectories looking for a file with the specified name
        string path = findFile(_applicationPath + _extdatpaths[i], name);
        if (!path.empty()) return System::canonicalPath(path);
    }

    // if we reach here, the resource wasn't found
    throw FATALERROR("Could not locate external resource '" + name + "'"
                     "\nDownload external resources from www.skirt.ugent.be using downloadextdat.sh");
}

////////////////////////////////////////////////////////////////////
