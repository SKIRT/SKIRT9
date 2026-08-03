/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "StoredColumns.hpp"
#include "FatalError.hpp"
#include "StringUtils.hpp"
#include "System.hpp"
#include <cstring>

////////////////////////////////////////////////////////////////////

namespace
{
    // the alternate interpretations for 8-byte items in the stored columns format
    union ScolItem
    {
        double doubleType;
        size_t sizeType;
        char stringType[8];
    };
    const size_t itemSize = sizeof(ScolItem);

    static_assert((sizeof(size_t) == 8) & (sizeof(double) == 8) & (itemSize == 8),
                  "Cannot properly declare union for items in stored columns format");
}

////////////////////////////////////////////////////////////////////

void StoredColumns::open(string filepath)
{
    // verify that we don't have an associated map
    if (!_filePath.empty()) throw FATALERROR("Another column file is already associated with this instance");

    // remember the path for the input file
    _filePath = filepath;

    // acquire a memory map for the file; the function returns zeros if the memory map cannot be created
    auto map = System::acquireMemoryMap(_filePath);
    if (!map.first) throw FATALERROR("Cannot acquire memory map for file: " + _filePath);
    const ScolItem* currentItem = static_cast<const ScolItem*>(map.first);

    // verify the name tag, the Endianness tag, and the extra zero
    if (memcmp("SKIRT X\n", currentItem++->stringType, itemSize) || currentItem++->sizeType != 0x010203040A0BFEFF
        || currentItem++->sizeType != 0)
        throw FATALERROR("File does not have stored columns file format: " + _filePath);

    // get the number of rows and columns
    size_t numRows = currentItem++->sizeType;
    _numColumns = currentItem++->sizeType;

    // get the column names and unit strings
    for (size_t i = 0; i != _numColumns; ++i)
        _columnNames.push_back(StringUtils::squeeze(string(currentItem++->stringType, itemSize)));
    for (size_t i = 0; i != _numColumns; ++i)
        _columnUnits.push_back(StringUtils::squeeze(string(currentItem++->stringType, itemSize)));

    // initialize the row pointers
    _nextRow = &currentItem->doubleType;
    currentItem += _numColumns * numRows;
    _endRow = &currentItem->doubleType;

    // verify the end-of-file tag
    if (memcmp("SCOLEND\n", currentItem->stringType, itemSize))
        throw FATALERROR("End-of-file does not match expected number of values: " + _filePath);
}

////////////////////////////////////////////////////////////////////

void StoredColumns::close()
{
    if (!_filePath.empty()) System::releaseMemoryMap(_filePath);
    _filePath.clear();
    _columnNames.clear();
    _columnUnits.clear();
    _numColumns = 0;
    _nextRow = nullptr;
    _endRow = nullptr;
}

////////////////////////////////////////////////////////////////////
