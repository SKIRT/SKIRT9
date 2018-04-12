/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef BUILDINFO_HPP
#define BUILDINFO_HPP

#include "Basics.hpp"

////////////////////////////////////////////////////////////////////

/** This class provides information about the current build, such as the build time. This
    information is provided to the source code at build time by CMake. */
class BuildInfo final
{
public:
    /** Returns the time of current build as a string formatted for human consumption. */
    static string timestamp();

    /** Returns the version string configured by CMake. */
    static string projectVersion();

    /** Returns a description of the code version from which this executable was built,
        e.g. the git commit hash. */
    static string codeVersion();
};

////////////////////////////////////////////////////////////////////

#endif
