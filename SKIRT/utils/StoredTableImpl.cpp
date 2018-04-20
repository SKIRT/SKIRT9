/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "StoredTableImpl.hpp"
#include "FatalError.hpp"
#include "StaticPaths.hpp"
#include "StringUtils.hpp"
#include "System.hpp"

////////////////////////////////////////////////////////////////////

namespace
{
    // the alternate interpretations for 8-byte items in the stored table format
    union Item
    {
        double doubleType;
        size_t sizeType;
        char stringType[8];
    };
    const size_t itemSize = sizeof(Item);

    static_assert((sizeof(size_t) == 8) & (sizeof(double) == 8) & (itemSize == 8),
                  "Cannot properly declare union for items in stored table format");
}

////////////////////////////////////////////////////////////////////

void StoredTable_Impl::open(size_t numAxes, string filename, string quantity,
                            string& filePath,
                            const double** axBeg, const double** qtyBeg,
                            size_t* axLen, size_t* qtyStep,
                            bool* axLog, bool* qtyLog)
{
    // add the mandatory filename extension if needed
    if (!StringUtils::endsWith(filename, ".stab")) filename += ".stab";

    // retrieve the canonical path for the resource; the function throws a fatal error if the resource cannot be found
    filePath = StaticPaths::resource(filename);

    // acquire a memory map for the resource; the function returns zeros if the memory map cannot be created
    auto map = System::acquireMemoryMap(filePath);
    if (!map.first) throw FATALERROR("Cannot acquire memory map for resource file: " + filePath);
    const Item* currentItem = static_cast<const Item*>(map.first);

    // verify the name tag and the Endianness tag
    if (memcmp("SKIRT X\n", currentItem++->stringType, itemSize) || currentItem++->sizeType != 0x010203040A0BFEFF)
        throw FATALERROR("Resource file does not have stored table format: " + filePath);

    // verify the number of axes
    if (currentItem++->sizeType != numAxes)
        throw FATALERROR("Number of axes in stored table file does not match: " + filePath);

    // skip the axes names and units, and store the interpolation scale
    for (size_t i = 0; i<numAxes; ++i) currentItem++;
    for (size_t i = 0; i<numAxes; ++i) currentItem++;
    for (size_t i = 0; i<numAxes; ++i) *axLog++ = (memcmp("log     ", currentItem++->stringType, itemSize) == 0);

    // store the grid length and a pointer to the first value for each axis
    for (size_t i = 0; i<numAxes; ++i)
    {
        size_t length = currentItem++->sizeType;
        *axLen++ = length;
        *axBeg++ = &currentItem->doubleType;
        currentItem += length;
    }

    // parse the requested quantity and unit from the specified string
    auto splitQuantity = StringUtils::split(quantity, ";");
    if (splitQuantity.size() != 2 || splitQuantity[0].size() < 1 || splitQuantity[1].size() < 1
                                  || splitQuantity[0].size() > 8 || splitQuantity[1].size() > 8)
        throw FATALERROR("Invalid tabulated quantity specification: " + quantity);
    string qntySpec = StringUtils::padRight(splitQuantity[0], itemSize);
    string unitSpec = StringUtils::padRight(splitQuantity[1], itemSize);

    // look for the appropriate quantity
    size_t numQnts = currentItem++->sizeType;
    size_t qntyIndex = numQnts;
    for (size_t i = 0; i<numQnts; ++i)
    {
        if (!memcmp(qntySpec.c_str(), currentItem++->stringType, itemSize)) qntyIndex = i;
    }
    if (qntyIndex==numQnts)
        throw FATALERROR("Tabulated quantity " + qntySpec + " is not in stored table " + filePath);

    // verify the corresponding unit
    for (size_t i = 0; i<numQnts; ++i)
    {
        if (i==qntyIndex)
        {
            if (memcmp(unitSpec.c_str(), currentItem++->stringType, itemSize))
                throw FATALERROR("Tabulated quantity does not have unit " + unitSpec + " in stored table " + filePath);
        }
        else currentItem++;
    }

    // store the corresponding interpolation scale
    for (size_t i = 0; i<numQnts; ++i)
    {
        if (i==qntyIndex) *qtyLog = (memcmp("log     ", currentItem++->stringType, itemSize) == 0);
        else currentItem++;
    }

    // calculate and store the pointer to the first quantity value, and store the number of quantities
    *qtyBeg = &currentItem->doubleType + qntyIndex;
    *qtyStep = numQnts;

}

////////////////////////////////////////////////////////////////////

void StoredTable_Impl::close(std::string filePath)
{
    if (!filePath.empty()) System::releaseMemoryMap(filePath);
}

////////////////////////////////////////////////////////////////////
