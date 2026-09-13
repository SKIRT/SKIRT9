/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef STOREDTABLEDICTIONARYIMPL_HPP
#define STOREDTABLEDICTIONARYIMPL_HPP

#include "Basics.hpp"
#include <unordered_map>
class SimulationItem;

////////////////////////////////////////////////////////////////////

/** This namespace contains private utilities for StoredTableDictionary. These non-template
    functions are seperated out into a seperate compilation unit (1) to avoid duplicating the code
    for each type of StoredTableDictionary and (2) to avoid propagating the dependencies (includes)
    of the implementation to all StoredTableDictionary users. */
namespace StoredTableDictionary_Impl
{
    /** Maps the name of a bundled stored table (including its mandatory ".stab" filename extension)
        to the byte offset, from the start of the bundle file, at which that stored table's data
        begins, and to that data's size in bytes. */
    using Index = std::unordered_map<string, std::pair<size_t, size_t>>;

    //============= Open and close =============

    /** This function performs the open() operation as described for the function with the same name
        in the StoredTableDictionary class template. It locates and memory-maps the specified bundle
        file (a tar archive in the POSIX ustar format, with the mandatory ".stabdict" filename
        extension added if needed), and populates \em index by walking the archive's headers,
        recording the byte offset and size of every regular-file member. Members whose base name
        starts with "._" (macOS AppleDouble sidecar files, which some tools add automatically when
        creating an archive) are silently skipped, as are non-regular-file entries such as directory
        or extended-header entries. The function throws a fatal error if the resource cannot be
        found, if a memory map cannot be acquired for it, or if it does not contain any usable
        members. */
    void open(const SimulationItem* item, string filename, string& filePath, Index& index);

    /** This function performs the close() operation as described for the destructor of the
        StoredTableDictionary class template. It receives the canonical path to the associated bundle
        file, or the empty string if no association exists. */
    void close(const string& filePath);

    //============= Query =============

    /** This function returns true if \em index has an entry for \em name (with the mandatory
        ".stab" filename extension added if needed), or false if not. */
    bool has(const Index& index, string name);

    /** This function returns the byte offset recorded in \em index for \em name (with the mandatory
        ".stab" filename extension added if needed). The \em filePath argument is used only to
        construct the error message thrown when \em name is not found in \em index. */
    size_t locate(const Index& index, const string& filePath, string name);
}

////////////////////////////////////////////////////////////////////

#endif
