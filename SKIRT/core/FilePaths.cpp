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

    // the expected resource pack names
    vector<string> _expectedPacks;

    // the expected resource pack version numbers:  <packname, version>
    std::unordered_map<string, int> _expectedPackVersions;

    // the installed resource pack version numbers:  <packname, version>
    std::unordered_map<string, int> _installedPackVersions;

    // relative paths to check for presence of built-in resources
    const char* _intpaths[] = {"../../../git/SKIRT/resources", "../../../../git/SKIRT/resources"};
    const int _Nintpaths = sizeof(_intpaths) / sizeof(const char*);

    // relative paths to check for presence of external resources
    const char* _extpaths[] = {"../../../resources", "../../../../resources"};
    const int _Nextpaths = sizeof(_extpaths) / sizeof(const char*);

    // recursively searches the given directory and its subdirectories for any files,
    // adds the corresponding canonical paths to the resource dictionary,
    // and processes any "version.txt" files encountered
    void findResourcesIn(string directory)
    {
        // search the current level
        for (const string& filename : System::filesInDirectory(directory))
        {
            if (!_resourcePaths.count(filename))
                _resourcePaths.emplace(filename, System::canonicalPath(StringUtils::joinPaths(directory, filename)));

            // remember the name and version number for each installed resource pack
            if (filename == "version.txt")
            {
                auto segments = StringUtils::split(directory, "_");
                if (segments.size() > 1)
                {
                    string packname = segments[segments.size() - 1];
                    int packversion = 0;
                    try
                    {
                        std::ifstream versionfile(StringUtils::joinPaths(directory, filename));
                        versionfile >> packversion;
                    }
                    catch (...)
                    {};

                    if (packversion > 0 && !_installedPackVersions.count(packname))
                        _installedPackVersions.emplace(packname, packversion);
                }
            }
        }

        // search the subdirectories
        for (const string& subdirname : System::dirsInDirectory(directory))
        {
            findResourcesIn(StringUtils::joinPaths(directory, subdirname));
        }
    }

    // parses the list of expected resource packs, if present, and stores the results
    // !! assumes that the dictionary with resource paths has already been populated !!
    void readExpectedPacks()
    {
        if (_resourcePaths.count("ExpectedResources.txt"))
        {
            try
            {
                std::ifstream expectedfile(_resourcePaths.at("ExpectedResources.txt"));
                while (expectedfile.good())
                {
                    string packname;
                    int packversion = 0;
                    expectedfile >> packname >> packversion;
                    if (!packname.empty() && packversion > 0 && !_expectedPackVersions.count(packname))
                    {
                        _expectedPacks.emplace_back(packname);
                        _expectedPackVersions.emplace(packname, packversion);
                    }
                }
            }
            catch (...)
            {};
        }
    }

    // populates the dictionary with resource paths and gathers resource pack information,
    // or throws a fatal error if there is a problem
    void findResources()
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
            if (System::isDir(test)) findResourcesIn(test);
        }

        // at this point we should have found at least some internal resources
        if (_resourcePaths.empty())
            throw FATALERROR("Could not locate built-in resources relative to '" + executableDirPath + "'");

        // iterate over the relative paths for external resources
        for (int i = 0; i < _Nextpaths; i++)
        {
            string test = StringUtils::joinPaths(executableDirPath, _extpaths[i]);
            if (System::isDir(test)) findResourcesIn(test);
        }

        // parse the list of expected resource packs, if present
        readExpectedPacks();
    }

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

FilePaths::FilePaths(SimulationItem* parent)
{
    parent->addChild(this);
}

////////////////////////////////////////////////////////////////////

void FilePaths::setupSelfBefore()
{
    SimulationItem::setupSelfBefore();
    std::call_once(_initialized, findResources);
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

string FilePaths::input(string name) const
{
    if (StringUtils::isAbsolutePath(name)) return name;
    return _inputPath + name;
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

string FilePaths::output(string name) const
{
    return _outputPath + _outputPrefix + "_" + name;
}

////////////////////////////////////////////////////////////////////

string FilePaths::resource(string name)
{
    std::call_once(_initialized, findResources);

    if (!_resourcePaths.count(name))
    {
        throw FATALERROR("Could not locate resource '" + name
                         + "'"
                           "\nRun ./downloadResources.sh in the ~/SKIRT/git directory"
                           "\nOr download additional resources from www.skirt.ugent.be");
    }
    return _resourcePaths.at(name);
}

////////////////////////////////////////////////////////////////////

string FilePaths::resourceName(string type, const vector<string>& segments)
{
    std::call_once(_initialized, findResources);

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

vector<string> FilePaths::expectedPacks()
{
    std::call_once(_initialized, findResources);

    return _expectedPacks;
}

////////////////////////////////////////////////////////////////////

int FilePaths::expectedPackVersion(string name)
{
    std::call_once(_initialized, findResources);

    return _expectedPackVersions.count(name) ? _expectedPackVersions.at(name) : 0;
}

////////////////////////////////////////////////////////////////////

int FilePaths::installedPackVersion(string name)
{
    std::call_once(_initialized, findResources);

    return _installedPackVersions.count(name) ? _installedPackVersions.at(name) : 0;
}

////////////////////////////////////////////////////////////////////
