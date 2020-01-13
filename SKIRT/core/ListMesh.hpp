/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef LISTMESH_HPP
#define LISTMESH_HPP

#include "TabulatedMesh.hpp"

////////////////////////////////////////////////////////////////////

/** ListMesh is a subclass of the MoveableMesh class. It represents a one-dimensional mesh over the
    unit interval [0,1] with mesh points that are fully specified inside the configuration file
    (i.e. without referring to an input file). It is intended for use in cases where there are just
    a mesh border points, but nothing keeps the user from specifying a long list. The \em numBins
    property of the Mesh base class is overridden to match the number of bins defined by the user.

    The configured list specifies the mesh points (i.e. the border points separating the mesh bins)
    in arbitrary order and in arbitary units. If the smallest point is not zero, an extra zero mesh
    point is inserted. The mesh is scaled so that the largest point in the list is mapped to unity.
    */
class ListMesh : public TabulatedMesh
{
    ITEM_CONCRETE(ListMesh, TabulatedMesh, "a mesh specified inside the configuration file")

        PROPERTY_DOUBLE_LIST(points, "the mesh border points")
        ATTRIBUTE_MIN_VALUE(points, "[0")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

    /** This function loads the mesh border points from the configured file and returns them, in
        arbitrary order and with arbitrary scaling. */
    vector<double> getMeshBorderPoints() const override;
};

////////////////////////////////////////////////////////////////////

#endif
