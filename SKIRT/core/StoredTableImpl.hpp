/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef STOREDTABLEIMPL_HPP
#define STOREDTABLEIMPL_HPP

#include "Array.hpp"
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

    /** This function performs the close() operation as described for the destructor of the
        StoredTable class template. It receives the canonical path to the associated resource file,
        or the empty string if no association exists. */
    void close(string filePath);
}

////////////////////////////////////////////////////////////////////

#endif
