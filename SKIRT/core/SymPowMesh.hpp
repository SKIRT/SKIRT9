/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef SYMPOWMESH_HPP
#define SYMPOWMESH_HPP

#include "MoveableMesh.hpp"

//////////////////////////////////////////////////////////////////////

/** The SymPowMesh class represents meshes on the unit interval \f$[0,1]\f$ with a symmetric
    power-law distribution of the mesh points. This distribution is such that the bin sizes form a
    geometric series, starting from the innermost bin and moving outwards symmetrically. If the
    number of bins is odd, there is one innermost bin; if it is even, there are two equal-size
    innermost bins. The actual distribution is characterized by a single parameter: the bin width
    ratio between the outermost and the innermost bins. This ratio can be larger than one (in which
    case the bin size increases when moving outwards) or smaller than one (in which case the bin
    size decrease from the centre to the edge of the interval). If the mesh has only one bin, this
    bin spans the complete interval \f$[0,1]\f$. If the mesh has two bins, each bin spans exactly
    half of the interval \f$[0,1]\f$. */
class SymPowMesh : public MoveableMesh
{
    ITEM_CONCRETE(SymPowMesh, MoveableMesh, "a symmetric power-law mesh")

        PROPERTY_DOUBLE(ratio, "the bin width ratio between the outermost and the innermost bins")
        ATTRIBUTE_MIN_VALUE(ratio, "]0")
        ATTRIBUTE_DEFAULT_VALUE(ratio, "1")

    ITEM_END()

    //======================== Other Functions =======================

public:
    /** This function returns an array containing the mesh points. */
    Array mesh() const override;
};

//////////////////////////////////////////////////////////////////////

#endif
