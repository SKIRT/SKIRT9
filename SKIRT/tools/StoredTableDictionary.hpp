/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef STOREDTABLEDICTIONARY_HPP
#define STOREDTABLEDICTIONARY_HPP

#include "StoredTable.hpp"
#include "StoredTableDictionaryImpl.hpp"

////////////////////////////////////////////////////////////////////

/** An instance of the StoredTableDictionary<N> class template provides access to a collection of
    stored tables (see StoredTable<N>) bundled into a single resource file, called a "dictionary",
    without requiring a separate memory map (and thus a separate open file) for each of the bundled
    tables.

    Motivation
    ----------

    A StoredTable<N> instance acquires a memory map on its associated stored table resource file for
    as long as the instance exists (see the StoredTable<N> class header for the rationale). This
    works well when the number of distinct stored table files involved stays modest. However, some
    use cases require selecting, at runtime, an a priori unknown but potentially large subset of
    tables out of a large family of similarly-shaped stored table resource files. Opening each of
    these as an indepedent StoredTable<N> would require a memory map, and thus an open file
    descriptor, per table, which can exceed operating system limits on the number of concurrently
    open (mapped) files.

    A dictionary file addresses this by bundling all of the individual stored table files into a
    single archive in the POSIX ustar format (i.e. the format written and read by the standard Unix
    "tar" tool), conventionally using the filename extension ".stabdict". Opening the dictionary
    acquires a single memory map on this bundle file and builds an in-memory index of the byte
    offset, within the bundle, of each of its members. A StoredTable<N> can then be spawned for any
    of these members without any further file I/O: it shares the dictionary's own memory map (which
    is safe because the sharing is refcounted; see below) and simply lets the stored table data
    format parser start reading at the appropriate offset instead of at the start of a dedicated
    memory map.

    The StoredTableDictionary<N> class template
    --------------------------------------------

    The template parameter \em N has the same meaning as for StoredTable<N>. All stored tables
    bundled into a given dictionary must share the same number of axes and the same axes and
    quantity specification, both matching the values passed to the dictionary's constructor or its
    open() function (using the same syntax as the corresponding StoredTable<N> arguments); the
    dictionary itself does not verify this up front, but every stored table spawned from it does, so
    a dictionary bundling incompatible members would report the problem at the point where such a
    member is actually spawned rather than at the point where the dictionary itself is opened.

    The default constructor creates an invalid stored table dictionary instance. The alternate
    constructor and the open() function associate a particular dictionary resource file with the
    instance. The has() function tests, and the open() function taking a member name spawns, a
    StoredTable<N> for one of the bundled members, identified by the filename it would have if it
    were a stand-alone stored table resource file (the mandatory ".stab" filename extension is added
    if needed, just as for StoredTable<N>'s own open() function).

    Because the sharing of the dictionary's memory map with its spawned stored tables is refcounted
    at the operating-system-resource level (see System::acquireMemoryMap()), a StoredTable<N> spawned
    from a dictionary remains fully valid even after the dictionary instance itself has been
    destroyed, for as long as at least one reference to the underlying memory map remains open,
    whether held by the dictionary or by any table spawned from it (or, for that matter, by another
    dictionary or stored table happening to reference the very same file). There is thus no need to
    keep a dictionary instance around for longer than it takes to spawn the stored tables of
    interest, although doing so avoids the (mild) overhead of repeatedly re-acquiring the memory map
    if further stored tables are spawned from it later on. */
template<size_t N> class StoredTableDictionary
{
public:
    /** The default constructor constructs an invalid stored table dictionary instance. The user
        must call the open() function (the one taking a dictionary filename) to associate the
        instance with a particular dictionary resource file. Calling any of the other functions
        before doing so results in undefined behavior (usually a crash). */
    StoredTableDictionary() {}

    /** This alternate constructor constructs a stored table dictionary instance and immediately
        associates a given dictionary resource file with it by calling the open() function taking a
        dictionary filename. Refer to that function for a description of the arguments and of its
        operation. */
    StoredTableDictionary(const SimulationItem* item, string filename, string axes, string quantity)
    {
        open(item, filename, axes, quantity);
    }

    /** The copy and move constructors and assignment operators are deleted; a dictionary is meant to
        be constructed once and then used in place (for example as a member of a longer-lived
        registry object). */
    StoredTableDictionary(const StoredTableDictionary&) = delete;
    StoredTableDictionary& operator=(const StoredTableDictionary&) = delete;
    StoredTableDictionary(StoredTableDictionary&&) = delete;
    StoredTableDictionary& operator=(StoredTableDictionary&&) = delete;

    /** This function associates a given dictionary resource file with the dictionary instance. If
        such an association already exists, the behavior is undefined. The \em item argument
        specifies a simulation item in the hierarchy of the caller (usually the caller itself), used
        to retrieve an appropriate logger for the stored tables later spawned from this dictionary.
        The \em filename argument specifies the filename of the resource, without any directory
        segments; it must have the ".stabdict" filename extension, which is added if needed.

        The \em axes and \em quantity arguments specify the axes and quantity specification shared
        by every stored table bundled into this dictionary, using the same syntax as the
        corresponding arguments of StoredTable<N>::open(); they are recorded for later use by open()
        (the one taking a member name) and are not verified against the dictionary's contents at this
        point.

        This function (1) locates the dictionary resource file, (2) acquires a memory map on it, and
        (3) builds an in-memory index of its members. If any of these steps fail, or if the
        dictionary turns out to have no usable members, the function throws a fatal error. */
    void open(const SimulationItem* item, string filename, string axes, string quantity)
    {
        StoredTableDictionary_Impl::open(item, filename, _filePath, _index);
        _axes = axes;
        _quantity = quantity;
    }

    /** The destructor releases this instance's reference on the dictionary's memory map, if there is
        one. As described in the class header, any stored table previously spawned from this
        dictionary remains valid regardless. */
    ~StoredTableDictionary() { StoredTableDictionary_Impl::close(_filePath); }

    // ================== Querying and spawning ==================

public:
    /** This function returns true if this dictionary has a member with the given name, or false if
        not. The \em name argument specifies the filename the member would have if it were a
        stand-alone stored table resource file, without any directory segments; the mandatory ".stab"
        filename extension is added if needed. */
    bool has(string name) const { return StoredTableDictionary_Impl::has(_index, name); }

    /** This function spawns and returns, by value, a StoredTable<N> for the member of this
        dictionary with the given name, sharing this dictionary's memory map rather than mapping a
        dedicated file (see the class header for the consequences). If this dictionary has no member
        with the given name, this function throws a fatal error. Unlike StoredTable<N>::open(), this
        function does not log anything on success, since the dictionary itself already logged a
        single message when it was opened.

        The \em name argument specifies the filename the member would have if it were a stand-alone
        stored table resource file, without any directory segments; the mandatory ".stab" filename
        extension is added if needed.

        The number of axes in the spawned stored table must match the template parameter \em N, and
        its axes and quantity must match the specifications passed to this dictionary's own
        constructor or open() function, exactly as described for StoredTable<N>::open(). The optional
        \em clampFirstAxis argument has the same meaning as for StoredTable<N>::open(). */
    StoredTable<N> open(string name, bool clampFirstAxis = true) const
    {
        size_t byteOffset = StoredTableDictionary_Impl::locate(_index, _filePath, name);
        return StoredTable<N>(_filePath, byteOffset, _filePath + ":" + name, _axes, _quantity, clampFirstAxis);
    }

    // ================== Data members ==================

private:
    string _filePath;                          // the canonical path to the associated dictionary file
    string _axes;                              // the axes specification shared by every bundled table
    string _quantity;                          // the quantity specification shared by every bundled table
    StoredTableDictionary_Impl::Index _index;  // maps a member name to its byte offset and size
};

////////////////////////////////////////////////////////////////////

#endif
