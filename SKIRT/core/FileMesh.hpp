/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef FILEMESH_HPP
#define FILEMESH_HPP

#include "TabulatedMesh.hpp"

////////////////////////////////////////////////////////////////////

/** FileMesh is a subclass of the MoveableMesh class. It represents a one-dimensional mesh over the
    unit interval [0,1] with mesh points that are read from a file. The \em numBins property of the
    Mesh base class is overridden to match the number of bins defined by the file.

    The input text file contains the mesh points (i.e. the border points separating the mesh bins),
    one point per line, in arbitrary order. The points may be given in arbitary units. If the
    smallest point is not zero, an extra zero mesh point is inserted. The mesh is scaled so that
    the largest point in the file is mapped to unity. */
class FileMesh : public TabulatedMesh
{
    ITEM_CONCRETE(FileMesh, TabulatedMesh, "a mesh read from a file")

        PROPERTY_STRING(filename, "the name of the file with the mesh border points")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

    /** This function loads the mesh border points from the configured file and returns them, in
        arbitrary order and with arbitrary scaling. */
    vector<double> getMeshBorderPoints() const override;
};

////////////////////////////////////////////////////////////////////

#endif
