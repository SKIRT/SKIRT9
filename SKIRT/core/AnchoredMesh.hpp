/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef ANCHOREDMESH_HPP
#define ANCHOREDMESH_HPP

#include "Mesh.hpp"

//////////////////////////////////////////////////////////////////////

/** The AnchoredMesh class represents meshes that are anchored to the origin; i.e. they can be
    scaled but not moved. This class does not add any functionality, it just serves to
    differentiate anchored meshes (which can be meaningfully employed only in the radial direction)
    and moveable meshes (which can be employed anywhere). */
class AnchoredMesh : public Mesh
{
    ITEM_ABSTRACT(AnchoredMesh, Mesh, "an anchored mesh")
    ITEM_END()
};

//////////////////////////////////////////////////////////////////////

#endif
