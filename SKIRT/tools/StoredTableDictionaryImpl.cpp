/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "StoredTableDictionaryImpl.hpp"
#include "FatalError.hpp"
#include "FilePaths.hpp"
#include "Log.hpp"
#include "StringUtils.hpp"
#include "System.hpp"
#include <cstring>

////////////////////////////////////////////////////////////////////

namespace
{
    // a POSIX ustar header block, and the associated data blocks, are always a multiple of this size
    constexpr size_t blockSize = 512;

    // byte offsets and sizes, within a ustar header block, of the fields we need
    constexpr size_t nameOffset = 0, nameSize = 100;
    constexpr size_t sizeOffset = 124, sizeSize = 12;
    constexpr size_t typeflagOffset = 156;
    constexpr size_t prefixOffset = 345, prefixSize = 155;

    // returns the given fixed-size, NUL-padded (but not necessarily NUL-terminated) header field as
    // a string truncated at the first NUL character, or spanning the full field if there is none
    string fieldToString(const char* field, size_t fieldSize)
    {
        return string(field, std::find(field, field + fieldSize, '\0'));
    }

    // parses the given fixed-size header field as an ASCII-octal, NUL- or space-terminated number
    size_t parseOctalField(const char* field, size_t fieldSize)
    {
        return static_cast<size_t>(std::strtoull(string(field, fieldSize).c_str(), nullptr, 8));
    }

    // adds the mandatory filename extension to the given name if needed
    string normalize(string name)
    {
        if (!StringUtils::endsWith(name, ".stab")) name += ".stab";
        return name;
    }
}

////////////////////////////////////////////////////////////////////

void StoredTableDictionary_Impl::open(const SimulationItem* item, string filename, string& filePath, Index& index)
{
    // add the mandatory filename extension if needed
    if (!StringUtils::endsWith(filename, ".stabdict")) filename += ".stabdict";

    // retrieve the full path for the resource; the function throws a fatal error if the resource cannot be found
    filePath = FilePaths::resource(filename);

    // acquire a memory map for the file; the function returns zeros if the memory map cannot be created
    auto map = System::acquireMemoryMap(filePath);
    if (!map.first) throw FATALERROR("Cannot acquire memory map for file: " + filePath);
    const char* base = static_cast<const char*>(map.first);
    size_t length = map.second;

    // walk the sequence of ustar header blocks, each followed by that member's data rounded up to a
    // multiple of the block size; the archive ends at the first all-zero header block (or at the end
    // of the mapped region, if the terminating blocks are missing)
    static const char zeroBlock[blockSize] = {};
    size_t offset = 0;
    while (offset + blockSize <= length)
    {
        const char* header = base + offset;
        if (memcmp(header, zeroBlock, blockSize) == 0) break;

        string name = fieldToString(header + nameOffset, nameSize);
        string prefix = fieldToString(header + prefixOffset, prefixSize);
        if (!prefix.empty()) name = prefix + "/" + name;
        size_t size = parseOctalField(header + sizeOffset, sizeSize);
        char typeflag = header[typeflagOffset];
        size_t dataOffset = offset + blockSize;

        // keep only regular-file entries, skipping macOS AppleDouble sidecar files ("._name") that
        // some tools add automatically alongside every real member
        if (typeflag == '0' || typeflag == '\0')
        {
            auto slash = name.find_last_of('/');
            string basename = (slash == string::npos) ? name : name.substr(slash + 1);
            if (!StringUtils::startsWith(basename, "._")) index[name] = {dataOffset, size};
        }

        offset = dataOffset + (size + blockSize - 1) / blockSize * blockSize;
    }

    if (index.empty()) throw FATALERROR("No stored tables found in dictionary file: " + filePath);

    item->find<Log>()->info(item->type() + " opened stored table dictionary " + filePath + " ("
                            + std::to_string(index.size()) + " members)");
}

////////////////////////////////////////////////////////////////////

void StoredTableDictionary_Impl::close(const string& filePath)
{
    if (!filePath.empty()) System::releaseMemoryMap(filePath);
}

////////////////////////////////////////////////////////////////////

bool StoredTableDictionary_Impl::has(const Index& index, string name)
{
    return index.count(normalize(name)) != 0;
}

////////////////////////////////////////////////////////////////////

size_t StoredTableDictionary_Impl::locate(const Index& index, const string& filePath, string name)
{
    name = normalize(name);
    auto it = index.find(name);
    if (it == index.end()) throw FATALERROR("Stored table " + name + " not found in dictionary file: " + filePath);
    return it->second.first;
}

////////////////////////////////////////////////////////////////////
