/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "StoredTableImpl.hpp"
#include "FatalError.hpp"
#include "FilePaths.hpp"
#include "Log.hpp"
#include "StringUtils.hpp"
#include "System.hpp"
#include <cstring>

////////////////////////////////////////////////////////////////////

namespace
{
    // the alternate interpretations for 8-byte items in the stored table format
    union StabItem
    {
        double doubleType;
        size_t sizeType;
        char stringType[8];
    };
    const size_t itemSize = sizeof(StabItem);

    static_assert((sizeof(size_t) == 8) & (sizeof(double) == 8) & (itemSize == 8),
                  "Cannot properly declare union for items in stored table format");
}

////////////////////////////////////////////////////////////////////

namespace
{
    // returns true if the name part of the given specification string matches the given 8-byte item
    bool matchesName(string specification, const StabItem* nameItem)
    {
        // parse the name part from the specified string
        auto index = specification.find('(');
        if (index == string::npos || index == 0 || index > 8)
            throw FATALERROR("Invalid stored table name specification: " + specification);
        string nameString = StringUtils::padRight(specification.substr(0, index), itemSize);

        // compare with the item
        return memcmp(nameString.c_str(), nameItem->stringType, itemSize) == 0;
    }

    // returns true if the unit part of the given specification string matches the given 8-byte item
    bool matchesUnit(string specification, const StabItem* unitItem)
    {
        // parse the unit part from the specified string
        auto size = specification.size();
        auto index = specification.find('(');
        if (index == string::npos || size - index <= 2 || size - index > 10 || specification[size - 1] != ')')
            throw FATALERROR("Invalid stored table unit specification: " + specification);
        string unitString = StringUtils::padRight(specification.substr(index + 1, size - index - 2), itemSize);

        // compare with the item
        return memcmp(unitString.c_str(), unitItem->stringType, itemSize) == 0;
    }
}

////////////////////////////////////////////////////////////////////

void StoredTable_Impl::open(size_t numAxes, const SimulationItem* item, string filename, bool resource, string axes,
                            string quantity, string& filePath, const double** axBeg, const double** qtyBeg,
                            size_t* axLen, size_t* qtyStep, bool* axLog, bool* qtyLog)
{
    // add the mandatory filename extension if needed
    if (!StringUtils::endsWith(filename, ".stab")) filename += ".stab";

    if (resource)
    {
        // retrieve the full path for the resource; the function throws a fatal error if the resource cannot be found
        filePath = FilePaths::resource(filename);
    }
    else
    {
        // retrieve the full path for the input file
        filePath = item->find<FilePaths>()->input(filename);
    }

    // acquire a memory map for the file; the function returns zeros if the memory map cannot be created
    auto map = System::acquireMemoryMap(filePath);
    if (!map.first) throw FATALERROR("Cannot acquire memory map for file: " + filePath);
    const StabItem* currentItem = static_cast<const StabItem*>(map.first);

    // verify the name tag and the Endianness tag
    if (memcmp("SKIRT X\n", currentItem++->stringType, itemSize) || currentItem++->sizeType != 0x010203040A0BFEFF)
        throw FATALERROR("File does not have stored table format: " + filePath);

    // verify the number of axes
    if (currentItem++->sizeType != numAxes)
        throw FATALERROR("Number of axes in stored table file does not match: " + filePath);

    // split the axes specification into a list
    auto axesSpecs = StringUtils::split(axes, ",");
    if (axesSpecs.size() != numAxes)
        throw FATALERROR("Number of axes in stored table axes specification does not match: " + axes);

    // verify the axes names
    for (size_t i = 0; i < numAxes; ++i)
    {
        if (!matchesName(axesSpecs[i], currentItem++))
            throw FATALERROR("Axis " + std::to_string(i) + " does not have name " + axesSpecs[i] + " in stored table "
                             + filePath);
    }

    // verify the axes units
    for (size_t i = 0; i < numAxes; ++i)
    {
        if (!matchesUnit(axesSpecs[i], currentItem++))
            throw FATALERROR("Axis " + std::to_string(i) + " does not have unit " + axesSpecs[i] + " in stored table "
                             + filePath);
    }

    // store the interpolation scale for each axis
    for (size_t i = 0; i < numAxes; ++i)
    {
        *axLog++ = (memcmp("log     ", currentItem++->stringType, itemSize) == 0);
    }

    // store the grid length and a pointer to the first value for each axis
    for (size_t i = 0; i < numAxes; ++i)
    {
        size_t length = currentItem++->sizeType;
        *axLen++ = length;
        *axBeg++ = &currentItem->doubleType;
        currentItem += length;
    }

    // look for the appropriate quantity
    size_t numQties = currentItem++->sizeType;
    size_t qtyIndex = numQties;
    for (size_t i = 0; i < numQties; ++i)
    {
        if (matchesName(quantity, currentItem++)) qtyIndex = i;
    }
    if (qtyIndex == numQties)
        throw FATALERROR("Tabulated quantity " + quantity + " is not in stored table " + filePath);

    // verify the corresponding unit
    for (size_t i = 0; i < numQties; ++i)
    {
        if (i == qtyIndex)
        {
            if (!matchesUnit(quantity, currentItem++))
                throw FATALERROR("Tabulated quantity does not have unit " + quantity + " in stored table " + filePath);
        }
        else
            currentItem++;
    }

    // store the corresponding interpolation scale
    for (size_t i = 0; i < numQties; ++i)
    {
        if (i == qtyIndex)
            *qtyLog = (memcmp("log     ", currentItem++->stringType, itemSize) == 0);
        else
            currentItem++;
    }

    // calculate and store the pointer to the first quantity value, and store the number of quantities
    *qtyBeg = &currentItem->doubleType + qtyIndex;
    *qtyStep = numQties;

    // log success, unless the same thread already successfully opened the same file for the same item
    thread_local const SimulationItem* previousItem = nullptr;
    thread_local string previousFilePath;
    if (item != previousItem || filePath != previousFilePath)
    {
        item->find<Log>()->info(item->type() + " opened stored table " + filePath);
        previousItem = item;
        previousFilePath = filePath;
    }
}

////////////////////////////////////////////////////////////////////

void StoredTable_Impl::close(std::string filePath)
{
    if (!filePath.empty()) System::releaseMemoryMap(filePath);
}

////////////////////////////////////////////////////////////////////
