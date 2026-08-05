/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef SYMLOGMESH_HPP
#define SYMLOGMESH_HPP

#include "Mesh.hpp"

//////////////////////////////////////////////////////////////////////

/** The SymLogMesh class represents a symmetric, centrally-anchored, logarithmic mesh. The mesh is
    actually composed of two centrally-anchored logarithmic meshes mirrored around the central
    point. The user-configurable parameter \f$t_\text{c}\f$ specifies the width of the central
    bin(s), and the remaining bins are distributed logarithmically over the intervals left and
    right of the central bin(s).

    If the total number of bins is odd, the single central bin covers the interval
    \f$[0.5-t_\text{c}/2, 0.5+t_\text{c}/2]\f$. If the number of bins is even, there are two
    central bins respectively covering the intervals \f$[0.5-t_\text{c}/2, 0.5]\f$ and \f$[0.5,
    0.5+t_\text{c}/2]\f$. The remaining bins are distributed symmetrically over the intervals
    \f$[0, 0.5-t_\text{c}/2]\f$ and \f$[0.5+t_\text{c}/2, 1]\f$, left and right of the central
    bin(s), using logarithmic spacing, with the smaller bins near the center.

    If the mesh has only one bin, the value of \f$t_\text{c}\f$ is ignored and the single bin spans
    the complete interval \f$[0,1]\f$. Similarly, if the mesh has only two bins, they respectively
    cover \f$[0,0.5]\f$ and \f$[0.5,1]\f$. */
class SymLogMesh : public Mesh
{
    ITEM_CONCRETE(SymLogMesh, Mesh, "a symmetric logarithmic mesh")

        PROPERTY_DOUBLE(centralBinFraction, "the central bin width fraction")
        ATTRIBUTE_MIN_VALUE(centralBinFraction, "]0")
        ATTRIBUTE_MAX_VALUE(centralBinFraction, "1[")
        ATTRIBUTE_DEFAULT_VALUE(centralBinFraction, "1e-3")

    ITEM_END()

    //======================== Other Functions =======================

public:
    /** This function returns an array containing the mesh points. */
    Array mesh() const override;
};

//////////////////////////////////////////////////////////////////////

#endif
