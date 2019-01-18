/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef FILEPATHS_HPP
#define FILEPATHS_HPP

#include "SimulationItem.hpp"

////////////////////////////////////////////////////////////////////

/** The FilePaths class manages the paths for the input and output files of a simulation, and for
    the resources included with the code or provided externally. */
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

    //======== Setters & Getters for Discoverable Attributes =======

public:
    /** Sets the (absolute or relative) path for input files. An empty string (the default value)
        means the current directory. */
    void setInputPath(string value);

    /** Returns the (absolute or relative) path for input files. */
    string inputPath() const;

    /** Sets the (absolute or relative) path for output files. An empty string (the default value)
        means the current directory. */
    void setOutputPath(string value);

    /** Returns the (absolute or relative) path for output files. */
    string outputPath() const;

    /** Sets the prefix for output file names; the default is empty (i.e. no prefix). */
    void setOutputPrefix(string value);

    /** Returns the prefix for output file names. */
    string outputPrefix() const;

    //======================== Other Functions =======================

public:
    /** This function returns the absolute canonical path for an input file with the specified
        name, relative to the input path returned by inputPath(). */
    string input(string name) const;

    /** This function returns the absolute canonical path for an output file with the specified
        name, relative to the output path returned by outputPath(). The prefix returned by
        outputPrefix() is inserted in front of the filename specified here. The prefix and the
        filename are separated by an underscore. */
    string output(string name) const;

    /** This function returns the absolute canonical path for a resource with the specified
        filename. The filename should \em not include any directory segments (just the base
        filename and filename extension). The function first searches resource files provided as
        part of the SKIRT build tree in the source code repository, and then searches resource
        files installed in a directory next to (and thus outside of) the SKIRT build tree. This
        mechanism allows to provide small and frequently-used resource files as part of the source
        code repository, while requiring larger resource files to be downloaded seperately.

        Specifically, the function searches the following directories, and all nested
        subdirectories inside these directories, recursively:

        - the \c resources directory inside the SKIRT build tree.
        - the \c resources directory (if any) next to the SKIRT \c git directory (i.e. outside of
          the build tree).

        The top-level directories are searched in the order listed above. The search order for the
        nested directories inside the top-level directories is unspecified. The first occurrence of
        the specified filename terminates the search. If the function cannot locate the specified
        resource, a fatal error is thrown. */
    static string resource(string name);

    /** This function returns the filename (without directory segments) for a resource with the
        specified type and with a filename including the specified segments. The list of resources
        scanned by this function is the same as that described for the resource() function.

        For a resource to be considered by this function, the end of its filename (including the
        filename extension) must match the specified type string, and the filename must also
        contain each of the specified segments. If no resources or multiple resources match these
        requirements, the function throws a fatal error. If a single resource matches, the function
        returns its filename. */
    static string resourceName(string type, const vector<string>& segments);

    //======================== Data Members ========================

private:
    string _inputPath;
    string _outputPath;
    string _outputPrefix;
};

////////////////////////////////////////////////////////////////////

#endif
