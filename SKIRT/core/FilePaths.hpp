/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef FILEPATHS_HPP
#define FILEPATHS_HPP

#include "SimulationItem.hpp"

////////////////////////////////////////////////////////////////////

/** The FilePaths class manages the paths for the input and output files of a simulation, and for
    the resources included with the code or provided externally.

    <b>Input and output files</b>

    A client program determines the input and output file paths and the prefix for output file
    names, presumably from the command line options. It then stores this information in the
    FilePaths instance held by the simulation, so that the simulation hierarchy can easily access
    it.

    <b>Resource files</b>

    The FilePaths class also offers static functions for locating the resource files used by the
    simulation hierarchy. These files can be provided as part of the build tree in the source code
    repository, or can be installed in a directory next to (and thus outside of) the build tree.
    This mechanism allows to provide small and frequently-used resource files as part of the source
    code repository, while requiring larger resource files to be downloaded seperately.

    Resource files are identified by their filename, without any directory or path information. In
    other words, a given resource file could be located in any of the supported directories (or
    nested subdirectory). On the other hand, this means there can be only a single resource with a
    particular name.

    Specifically, when first invoked, the FilePaths class builds a list of all available resource
    files by iterating over the following directories and all nested subdirectories inside these
    directories, recursively:

    - the \c resources directory inside the SKIRT build tree.

    - the \c resources directory (if any) next to the SKIRT \c git directory (i.e. outside of the
    build tree).

    The top-level directories are iterated in the order listed above. The iteration order for the
    nested subdirectories inside the top-level directories is unspecified. The first occurrence of
    a particular filename is stored; any subsequent occurrences of the same filename are ignored.

    <b>Resource packs</b>

    Finally, the FilePaths class offers static functions that support versioning for resource \em
    packs, which are usually installed in the resource directory outside of the build tree.

    A resource pack is a set of resource files contained in a resource subdirectory, called the
    pack directory, that contains a file named \c version.txt. The first token in this text file's
    contents is interpreted as the integer version number of the resource pack. The name of the
    resource pack is determined by the segment of the pack directory name after the last
    underscore. The resource files in the pack can reside in nested subdirectories.

    Furthermore, a resource file named \c ExpectedResources.txt can be provided as part of the
    build tree. This text file lists the names and version numbers of the resource packs expected
    to be installed, one pack name and corresponding version number per line.

    The FilePaths class offers functions to retrieve information on the expected and installed
    resource packs from these files. This information can be used by the client program to verify
    that the appropriate resource pack versions have been installed.

    */
class FilePaths : public SimulationItem
{
    //============= Construction - Setup - Destruction =============

public:
    /** This constructor creates a file path object that is hooked up as a child to the specified
        parent in the simulation hierarchy, so that it will automatically be deleted. The setup()
        function is \em not called by this constructor. */
    explicit FilePaths(SimulationItem* parent);

protected:
    /** This function determines and caches the resource file paths that can be returned by this
        class. This avoids repeated searches through the resource directories, and allows reporting
        any problems as early as possible in the program's lifecycle. */
    void setupSelfBefore() override;

    //======== Input and output file paths and names =======

public:
    /** Sets the (absolute or relative) path for input files. An empty string (the default value)
        means the current directory. */
    void setInputPath(string value);

    /** Returns the (absolute or relative) path for input files. */
    string inputPath() const;

    /** This function returns the absolute canonical path for an input file with the specified
        name, relative to the input path returned by inputPath(). */
    string input(string name) const;

    /** Sets the (absolute or relative) path for output files. An empty string (the default value)
        means the current directory. */
    void setOutputPath(string value);

    /** Returns the (absolute or relative) path for output files. */
    string outputPath() const;

    /** Sets the prefix for output file names; the default is empty (i.e. no prefix). */
    void setOutputPrefix(string value);

    /** Returns the prefix for output file names. */
    string outputPrefix() const;

    /** This function returns the absolute canonical path for an output file with the specified
        name, relative to the output path returned by outputPath(). The prefix returned by
        outputPrefix() is inserted in front of the filename specified here. The prefix and the
        filename are separated by an underscore. */
    string output(string name) const;

    //======================== Resource files =======================

public:
    /** This function returns the absolute canonical path for a resource file with the specified
        filename. The filename should \em not include any directory segments (just the base
        filename and filename extension). The function searches the list of available resource
        files as described in the class header. If the specified resource file cannot be located, a
        fatal error is thrown. */
    static string resource(string name);

    /** This function returns the filename (without directory segments) for a resource file with
        the specified type and with a filename including the specified segments. The function
        searches the list of available resource files as described in the class header.

        For a resource file to be considered by this function, the end of its filename (including
        the filename extension) must match the specified type string, and the filename must also
        contain each of the specified segments. If no or multiple resources files match these
        requirements, the function throws a fatal error. If a single resource file matches, the
        function returns its filename. */
    static string resourceName(string type, const vector<string>& segments);

    //======================== Resource packs =======================

public:
    /** This function returns a list of the names of all expected packs, in the order listed in the
        \c ExpectedResources.txt resource file, or the empty list if that resource file is missing
        or empty. */
    static vector<string> expectedPacks();

    /** This function returns the expected version number for the resource pack with the specified
        name, or zero if the name is not in the list of expected packs. */
    static int expectedPackVersion(string name);

    /** This function returns the version number for the installed resource pack with the specified
        name, or zero if the pack is not installed. */
    static int installedPackVersion(string name);

    //======================== Data Members ========================

private:
    string _inputPath;
    string _outputPath;
    string _outputPrefix;
};

////////////////////////////////////////////////////////////////////

#endif
