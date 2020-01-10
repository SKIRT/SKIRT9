/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef VORONOIMESHINTERFACE_HPP
#define VORONOIMESHINTERFACE_HPP

#include "Basics.hpp"
class VoronoiMeshSnapshot;

////////////////////////////////////////////////////////////////////

/** VoronoiMeshInterface is a pure interface. It provides access to the Voronoi mesh snapshot
    maintained by the object that implements the interface. */
class VoronoiMeshInterface
{
protected:
    /** The empty constructor for the interface. */
    VoronoiMeshInterface() {}

public:
    /** The empty destructor for the interface. */
    virtual ~VoronoiMeshInterface() {}

    /** This function must be implemented in a derived class. It returns a pointer to the Voronoi
        mesh snapshot maintained by the object that implements the interface. */
    virtual VoronoiMeshSnapshot* voronoiMesh() const = 0;
};

/////////////////////////////////////////////////////////////////////////////

#endif
