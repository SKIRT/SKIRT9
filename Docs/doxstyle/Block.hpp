/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef BLOCK_HPP
#define BLOCK_HPP

#include "Basics.hpp"

////////////////////////////////////////////////////////////////////

/** This class streamlines a single doxygen documentation block. */
class Block
{
public:
    /** Constructs a streamliner for the specified documentation block. */
    Block(const vector<string>& block);

    /** Constructs a streamliner for the documentation block specified as a portion of a larger
        source code chunk. */
    Block(const vector<string>& chunk, size_t first, size_t last);

    /** Returns the streamlined documentation block. */
    vector<string> streamlined();

private:
    /** The documentation block. */
    vector<string> _block;
};

////////////////////////////////////////////////////////////////////

#endif
