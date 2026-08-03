/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef LINMESH_HPP
#define LINMESH_HPP

#include "Mesh.hpp"

//////////////////////////////////////////////////////////////////////

/** The LinMesh class represents a mesh on the unit interval \f$[0,1]\f$ with a linear distribution
    of the mesh points. */
class LinMesh : public Mesh
{
    ITEM_CONCRETE(LinMesh, Mesh, "a linear mesh")
    ITEM_END()

    //======================== Other Functions =======================

public:
    /** This function returns an array containing the mesh points. */
    Array mesh() const override;
};

//////////////////////////////////////////////////////////////////////

#endif
