/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef TETRAMESHINTERFACE_HPP
#define TETRAMESHINTERFACE_HPP

#include "Basics.hpp"
class TetraMeshSnapshot;

////////////////////////////////////////////////////////////////////

/** TetraMeshInterface is a pure interface. It provides access to the Tetra mesh snapshot
    maintained by the object that implements the interface. */
class TetraMeshInterface
{
protected:
    /** The empty constructor for the interface. */
    TetraMeshInterface() {}

public:
    /** The empty destructor for the interface. */
    virtual ~TetraMeshInterface() {}

    /** This function must be implemented in a derived class. It returns a pointer to the Tetra
        mesh snapshot maintained by the object that implements the interface. */
    virtual TetraMeshSnapshot* tetraMesh() const = 0;
};

/////////////////////////////////////////////////////////////////////////////

#endif
