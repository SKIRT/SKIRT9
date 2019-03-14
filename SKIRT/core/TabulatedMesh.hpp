/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef TABULATEDMESH_HPP
#define TABULATEDMESH_HPP

#include "MoveableMesh.hpp"

////////////////////////////////////////////////////////////////////

/** TabulatedMesh is an abstract class for representing a one-dimensional mesh over the unit
    interval [0,1] with mesh points that are tabulated by the user. The \em numBins property of the
    Mesh base class is overridden to match the number of bins in the table defined by the user.

    The subclass must load the tabulated data, and this abstract class handles everything else.

    The loaded table contains the mesh points (i.e. the border points separating the mesh bins), in
    arbitrary order. The points may be given in arbitary units. If the smallest point is not zero,
    an extra zero mesh point is inserted. The mesh is scaled so that the largest point in the table
    is mapped to unity. */
class TabulatedMesh : public MoveableMesh
{
    ITEM_ABSTRACT(TabulatedMesh, MoveableMesh, "a mesh tabulated by the user")
        ATTRIBUTE_TYPE_DISPLAYED_IF(TabulatedMesh, "Level2")
    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function asks the subclass to load the mesh border points and precalculates
        the final mesh, sorting and scaling the values as described in the class header. */
    void setupSelfBefore() override;

    /** This function must be implemented in each subclass to return the mesh border points, in
        arbitrary order and with arbitrary scaling. */
    virtual vector<double> getMeshBorderPoints() const = 0;

    //======================== Other Functions =======================

public:
    /** This function returns an array containing the mesh points. */
    Array mesh() const override;

    //======================== Data Members ========================

private:
    Array _mesh;
};

////////////////////////////////////////////////////////////////////

#endif
