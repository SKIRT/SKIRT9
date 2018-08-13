/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef MOVEABLEMESH_HPP
#define MOVEABLEMESH_HPP

#include "Mesh.hpp"

//////////////////////////////////////////////////////////////////////

/** The MoveableMesh class represents meshes that can be moved (offset) in addition to being
    scaled. This class does not add any functionality, it just serves to differentiate anchored
    meshes (which can be meaningfully employed only in the radial direction) and moveable meshes
    (which can be employed anywhere). */
class MoveableMesh : public Mesh
{
    ITEM_ABSTRACT(MoveableMesh, Mesh, "a moveable mesh")
    ITEM_END()
};

//////////////////////////////////////////////////////////////////////

#endif
