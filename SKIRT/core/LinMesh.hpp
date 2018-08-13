/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef LINMESH_HPP
#define LINMESH_HPP

#include "MoveableMesh.hpp"

//////////////////////////////////////////////////////////////////////

/** The LinMesh class represents meshes on the unit interval \f$[0,1]\f$ with a linear distribution
    of the mesh points. */
class LinMesh : public MoveableMesh
{
    ITEM_CONCRETE(LinMesh, MoveableMesh, "a linear mesh")
    ITEM_END()

    //======================== Other Functions =======================

public:
    /** This function returns an array containing the mesh points. */
    Array mesh() const override;
};

//////////////////////////////////////////////////////////////////////

#endif
