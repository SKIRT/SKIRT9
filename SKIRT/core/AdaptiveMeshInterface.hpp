/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef ADAPTIVEMESHINTERFACE_HPP
#define ADAPTIVEMESHINTERFACE_HPP

#include "Basics.hpp"
class AdaptiveMeshSnapshot;

////////////////////////////////////////////////////////////////////

/** AdaptiveMeshInterface is a pure interface. It provides access to the adaptive mesh snapshot
    maintained by the object that implements the interface. */
class AdaptiveMeshInterface
{
protected:
    /** The empty constructor for the interface. */
    AdaptiveMeshInterface() {}

public:
    /** The empty destructor for the interface. */
    virtual ~AdaptiveMeshInterface() {}

    /** This function must be implemented in a derived class. It returns a pointer to the adaptive
        mesh snapshot maintained by the object that implements the interface. */
    virtual AdaptiveMeshSnapshot* adaptiveMesh() const = 0;
};

/////////////////////////////////////////////////////////////////////////////

#endif
