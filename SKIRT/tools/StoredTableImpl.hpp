/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef STOREDTABLEIMPL_HPP
#define STOREDTABLEIMPL_HPP

#include "Basics.hpp"
class SimulationItem;

////////////////////////////////////////////////////////////////////

/** This namespace contains private utilities for StoredTable. These non-template functions are
    seperated out into a seperate compilation unit (1) to avoid duplicating the code for each type
    of StoredTable and (2) to avoid propagating the dependencies (includes) of the implementation
    to all StoredTable users. */
namespace StoredTable_Impl
{
    //============= Open and close =============

    /** This function performs the open() operation as described for the function with the same name in the
        StoredTable class template. It receives references or pointers to all data members of the
        stored table instance, in addition to the input parameters of the open() function. */
    void open(size_t numAxes, const SimulationItem* item,   // input parameters
              string filename, bool resource,               //   "
              string axes, string quantity,                 //   "
              string& filePath,                             // output parameter by reference
              const double** axBeg, const double** qtyBeg,  // output parameters via pointers
              size_t* axLen, size_t* qtyStep,               //   "
              bool* axLog, bool* qtyLog);                   //   "

    /** This function performs the open() operation for a stored table whose data lives at a given
        byte offset within an already-resolved file, rather than being resolved from a bare resource
        or input filename as open() does. It is used by StoredTableDictionary to hand out a
        StoredTable that shares the dictionary's own memory-mapped bundle file instead of mapping a
        dedicated file of its own. The \em filePath argument must already be the canonical path of
        that bundle file (as previously obtained, for example, when the dictionary itself opened it);
        it is used to acquire and identify the shared memory map, so that this stored table shares
        and correctly refcounts that map with the dictionary and with every other stored table
        spawned from it. The \em byteOffset argument gives the position, in bytes from the start of
        that file, where this particular stored table's data begins. The \em label argument is used
        only to identify this particular table in error messages (typically combining the bundle
        file path and the member name); unlike open(), this function does not log anything on
        success, since the dictionary itself already logs a single message when it is opened. In all
        other respects, this function behaves as open() does: it verifies that the stored table
        matches all requirements, and stores relevant information in the given output parameters,
        throwing a fatal error if any of these steps fail. */
    void openAt(size_t numAxes,                               // input parameters
                string filePath, size_t byteOffset,           //   "
                string label,                                 //   "
                string axes, string quantity,                 //   "
                const double** axBeg, const double** qtyBeg,  // output parameters via pointers
                size_t* axLen, size_t* qtyStep,               //   "
                bool* axLog, bool* qtyLog);                   //   "

    /** This function performs the close() operation as described for the destructor of the
        StoredTable class template. It receives the canonical path to the associated resource file,
        or the empty string if no association exists. */
    void close(string filePath);
}

////////////////////////////////////////////////////////////////////

#endif
