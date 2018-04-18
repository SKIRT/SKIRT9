/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef FILEPATHS_HPP
#define FILEPATHS_HPP

#include "SimulationItem.hpp"

////////////////////////////////////////////////////////////////////

/** The FilePaths class manages the paths for the input and output files of a simulation. */
class FilePaths : public SimulationItem
{
    //============= Construction - Setup - Destruction =============

public:
    /** This constructor creates a file path object that is hooked up as a child to the specified
        parent in the simulation hierarchy, so that it will automatically be deleted. The setup()
        function is \em not called by this constructor. */
    explicit FilePaths(SimulationItem* parent);

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

    //======================== Data Members ========================

private:
    string _inputPath;
    string _outputPath;
    string _outputPrefix;
};

////////////////////////////////////////////////////////////////////

#endif
