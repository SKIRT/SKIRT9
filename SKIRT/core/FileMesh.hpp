/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef FILEMESH_HPP
#define FILEMESH_HPP

#include "MoveableMesh.hpp"

////////////////////////////////////////////////////////////////////

/** FileMesh is a subclass of the MoveableMesh class. It represents a one-dimensional mesh over the
    unit interval [0,1] with mesh points that are read from a file. The \em numBins property of the
    Mesh base class is overridden to match the number of bins defined by the file.

    The input text file contains the mesh points (i.e. the border points separating the mesh bins),
    one point per line, in arbitrary order. The points may be given in arbitary units. If the
    lowest point is not zero, an extra zero mesh point is inserted. The mesh is scaled so that the
    last point in the file is mapped to unity. */
class FileMesh : public MoveableMesh
{
    ITEM_CONCRETE(FileMesh, MoveableMesh, "a mesh read from a file")
        ATTRIBUTE_TYPE_DISPLAYED_IF(FileMesh, "Level2")

    PROPERTY_STRING(filename, "the name of the file with the mesh border points")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function reads the mesh points from the specified file. */
    void setupSelfBefore() override;

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
