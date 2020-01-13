/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "Chunk.hpp"
#include "Block.hpp"
#include "StringUtils.hpp"
#include <iostream>
#include <regex>

////////////////////////////////////////////////////////////////////

Chunk::Chunk() {}

////////////////////////////////////////////////////////////////////

void Chunk::readFromConsole()
{
    while (true)
    {
        string line;
        std::getline(std::cin, line);
        if (!std::cin) break;
        _chunk.push_back(line);
    }
}

////////////////////////////////////////////////////////////////////

void Chunk::writeToConsole()
{
    for (const string& line : _chunk)
    {
        std::cout << line << std::endl;
    }
}

////////////////////////////////////////////////////////////////////

namespace
{
    // Returns the index of the first string in the string list that matches the regular expression,
    // starting at the specified index, or the size of the list if none of the strings match
    size_t indexOf(const vector<string>& list, const std::regex& find, size_t index)
    {
        while (index < list.size())
        {
            if (std::regex_match(list[index], find)) break;
            index++;
        }
        return index;
    }
}

////////////////////////////////////////////////////////////////////

void Chunk::streamline()
{
    // some regular expressions
    std::regex startDox("\\s*/\\*\\*.*");  //   \s* (white space) /** (literal text) .* (anything)
    std::regex endDox(".*\\*/\\s*");       //    .* (anything) */ (literal text) \s* (white space)

    // loop over the lines of the source code
    size_t index = 0;
    while (index < _chunk.size())
    {
        // look for the start of a documentation block
        index = indexOf(_chunk, startDox, index);
        if (index >= _chunk.size()) break;

        // look for the end of the documentation block (which could be on the same line)
        size_t end = indexOf(_chunk, endDox, index);
        if (end >= _chunk.size()) break;

        // streamline the block
        Block styler(_chunk, index, end);
        vector<string> block = styler.streamlined();

        // replace the block in the source code
        _chunk.erase(_chunk.begin() + index, _chunk.begin() + end + 1);
        _chunk.insert(_chunk.begin() + index, block.begin(), block.end());
        index += block.size();

        // remove any empty lines following the block in the source code
        while (index < _chunk.size() && StringUtils::squeeze(_chunk[index]).empty())
            _chunk.erase(_chunk.begin() + index);
    }
}

////////////////////////////////////////////////////////////////////
