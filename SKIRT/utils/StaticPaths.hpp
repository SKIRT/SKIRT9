/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef STATICPATHS_HPP
#define STATICPATHS_HPP

#include "Basics.hpp"

////////////////////////////////////////////////////////////////////

/** The StaticPaths namespace retrieves global file paths, such as those for the resources included
    with the code or provided externally. */
namespace StaticPaths
{
    /** This function returns the absolute canonical path for a resource with the specified
        filename. The filename should \em not include any directory segments (just the base
        filename and filename extension). The function first looks for built-in resource files and
        then looks for externally provided resource files. This mechanism allows to provide small
        and frequently-used resource files as part of the SKIRT build tree in the source code
        repository, while requiring larger resource files to be downloaded seperately from the
        SKIRT web site using a shell script provided for this purpose.

        Specifically, the function searches the following directories, and all nested subdirectories
        inside these directories, recursively:

        - the \c resources directory inside the SKIRT build tree.
        - the \c resources directory (if any) next to the SKIRT \c git directory
          (i.e. outside of the build tree).

        The top-level directories are searched in the order listed above. The search order for the
        nested directories inside the top-level directories is unspecified. The first occurrence of
        the specified filename terminates the search. This means that one cannot replace a built-in
        resource by placing a file with the same name in an external directory.

        If the function cannot locate the specified resource, a fatal error is thrown. */
    string resource(string name);
}

////////////////////////////////////////////////////////////////////

#endif
