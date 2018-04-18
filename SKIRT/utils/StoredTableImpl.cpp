/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "StoredTableImpl.hpp"
#include "FatalError.hpp"
#include "System.hpp"

////////////////////////////////////////////////////////////////////

void StoredTable_Impl::open(size_t N, string filename, string quantity,
                            string& filePath,
                            const double** axBeg, const double** qtyBeg, size_t* axLen,
                            bool* axLog, bool* qtyLog)
{
    // TO DO...
    System::acquireMemoryMap(filePath);
}

////////////////////////////////////////////////////////////////////

void StoredTable_Impl::close(std::string filePath)
{
    if (!filePath.empty()) System::releaseMemoryMap(filePath);
}

////////////////////////////////////////////////////////////////////
