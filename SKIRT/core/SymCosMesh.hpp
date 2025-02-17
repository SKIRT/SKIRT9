/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef SYMCOSMESH_HPP
#define SYMCOSMESH_HPP

#include "Mesh.hpp"

//////////////////////////////////////////////////////////////////////

/** The SymCosMesh class represents a symmetric cosine mesh. The mesh points are spaced such that
    when the mesh is used to subdivide a spherical volume into conical slices by polar inclination
    angle \f$0 \le \theta \le \pi\f$, the resulting slices have equal volume. Specifically, the
    mesh is contructed as follows:

    - construct an equidistant linear mesh over the (inverted) interval \f$[1,-1]\f$ with the
    requested number of bins;

    - replace each mesh point by the arc-cosine of its value, resulting in a mesh with the intended
    spacing over the interval \f$[0,\pi]\f$;

    - divide each mesh point by \f$\pi\f$, rescaling the mesh to the unit interval \f$[0,1]\f$.

    The scaling in the latter step is required to conform to the standards for all Mesh classes.
    When the mesh is actually used to discretize an inclination angle, it will be scaled back to
    the range \f$[0,\pi]\f$. */
class SymCosMesh : public Mesh
{
    ITEM_CONCRETE(SymCosMesh, Mesh, "a symmetric cosine mesh")
    ITEM_END()

    //======================== Other Functions =======================

public:
    /** This function returns an array containing the mesh points. */
    Array mesh() const override;
};

//////////////////////////////////////////////////////////////////////

#endif
