/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "FilePaths.hpp"
#include "FatalError.hpp"
#include "StringUtils.hpp"
#include "System.hpp"
#include <mutex>
#include <unordered_map>

////////////////////////////////////////////////////////////////////

namespace
{
    // flag becomes true if the static paths have been initialized
    std::once_flag _initialized;

    // the resource paths:  <filename, complete_path>
    std::unordered_map<string, string> _resourcePaths;

    // relative paths to check for presence of built-in resources
    const char* _intpaths[] = {"../../../git/SKIRT/resources", "../../../../git/SKIRT/resources"};
    const int _Nintpaths = sizeof(_intpaths) / sizeof(const char*);

    // relative paths to check for presence of external resources
    const char* _extpaths[] = {"../../../resources", "../../../../resources"};
    const int _Nextpaths = sizeof(_extpaths) / sizeof(const char*);

    // recursively searches the given directory and its subdirectories for any files
    // and adds the corresponding canonical paths to the resource dictionary
    void findResourcePathsIn(string directory)
    {
        // search the current level
        for (const string& filename : System::filesInDirectory(directory))
        {
            if (!_resourcePaths.count(filename))
                _resourcePaths.emplace(filename, System::canonicalPath(StringUtils::joinPaths(directory, filename)));
        }

        // search the subdirectories
        for (const string& subdirname : System::dirsInDirectory(directory))
        {
            findResourcePathsIn(StringUtils::joinPaths(directory, subdirname));
        }
    }

    // populates the dictionary with resource paths, or throws an error if there is a problem
    void findResourcePaths()
    {
        // get the executable path (or the empty string in case of failure)
        string executableFilePath = System::executablePath();
        if (executableFilePath.empty()) throw FATALERROR("Could not determine path to executable");

        // get the location of the executable (i.e. the path to the containing directory)
        string executableDirPath = StringUtils::dirPath(executableFilePath);

        // iterate over the relative paths for built-in resources
        for (int i = 0; i < _Nintpaths; i++)
        {
            string test = StringUtils::joinPaths(executableDirPath, _intpaths[i]);
            if (System::isDir(test)) findResourcePathsIn(test);
        }

        // at this point we should have found at least some internal resources
        if (_resourcePaths.empty())
            throw FATALERROR("Could not locate built-in resources relative to '" + executableDirPath + "'");

        // iterate over the relative paths for external resources
        for (int i = 0; i < _Nextpaths; i++)
        {
            string test = StringUtils::joinPaths(executableDirPath, _extpaths[i]);
            if (System::isDir(test)) findResourcePathsIn(test);
        }
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
    std::call_once(_initialized, findResourcePaths);
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
    if (StringUtils::isAbsolutePath(name)) return name;
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
    // initialize the resource paths if needed
    std::call_once(_initialized, findResourcePaths);

    if (!_resourcePaths.count(name))
    {
        throw FATALERROR("Could not locate resource '" + name
                         + "'"
                           "\nDownload additional resources from www.skirt.ugent.be");
    }
    return _resourcePaths.at(name);
}

////////////////////////////////////////////////////////////////////

bool FilePaths::hasResource(std::string name)
{
    // initialize the resource paths if needed
    std::call_once(_initialized, findResourcePaths);

    return _resourcePaths.count(name) > 0;
}

////////////////////////////////////////////////////////////////////

namespace
{
    // returns true if the resource filename matches the requirements, false otherwise;
    // see documentation of the resourceName() function for more details
    bool matches(string resource, string type, const vector<string>& segments)
    {
        if (!StringUtils::endsWith(resource, type)) return false;
        for (string segment : segments)
        {
            if (!StringUtils::contains(resource, segment)) return false;
        }
        return true;
    }

    // returns a human-readable message describing the specified requirements;
    // see documentation of the resourceName() function for more details
    string message(string type, const vector<string>& segments)
    {
        string msg = "type '" + type + "'";
        if (!segments.empty())
        {
            msg += " with filename containing";
            for (string segment : segments) msg += " '" + segment + "',";
            msg.erase(msg.size() - 1, 1);
        }
        return msg;
    }
}

////////////////////////////////////////////////////////////////////

string FilePaths::resourceName(string type, const vector<string>& segments)
{
    // initialize the resource paths if needed
    std::call_once(_initialized, findResourcePaths);

    // look for matching resource filename
    string result;
    for (const auto& pair : _resourcePaths)
    {
        const string& resource = pair.first;
        if (matches(resource, type, segments))
        {
            // fail if there is ambiguity
            if (!result.empty()) throw FATALERROR("Multiple resources matching " + message(type, segments));
            result = resource;
        }
    }
    // fail if not found
    if (result.empty()) throw FATALERROR("No resources matching " + message(type, segments));
    return result;
}

////////////////////////////////////////////////////////////////////
