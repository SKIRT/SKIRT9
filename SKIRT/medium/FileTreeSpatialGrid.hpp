/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef FILETREESPATIALGRID_HPP
#define FILETREESPATIALGRID_HPP

#include "TreeSpatialGrid.hpp"

//////////////////////////////////////////////////////////////////////

/** FileTreeSpatialGrid is a concrete subclass of the TreeSpatialGrid class that constructs a
    hierarchical tree grid from the topological description loaded from a data file. In the most
    common use case, the data file has been created in a previous simulation through the
    TreeSpatialGridTopologyProbe. The tree grid in the current simulation will then be identical to
    the one in the previous simulation. This can be useful in situations where multiple simulations
    are being performed on input models with the same (or a very similar) spatial distribution of
    the transfer medium. Constructing a tree grid based on the medium density distribution takes
    time (because of the large number of density samples required), and the tree is likely to
    differ slightly between various runs (because the density is sampled at random positions, and
    the random number sequence varies because of parallelization). Loading the tree from the
    topology data file is much faster and guarantees that all simulations use identical tree grids.

    After a brief descriptive header, the input file contains lines with just a single integer
    number. The first line specifies the number of children for each nonleaf node (2 for a binary
    tree, 8 for an octtree, or 0 if the root node is not subdivided). The second line contains 1 if
    the root node is subdivided, or 0 if not. The following lines similarly contain 1 or 0
    indicating subdivision for any children of the preceding node, recursively, in a depth-first
    traversal of the tree.

    Note that the topology data stored in the input file is scale-free. When configuring the
    FileTreeSpatialGrid in the simulation loading the topology, the user must specify the extent of
    the spatial domain. */
class FileTreeSpatialGrid : public TreeSpatialGrid
{
    ITEM_CONCRETE(FileTreeSpatialGrid, TreeSpatialGrid, "a tree-based spatial grid loaded from a topology data file")
        ATTRIBUTE_TYPE_DISPLAYED_IF(FileTreeSpatialGrid, "Level2")

        PROPERTY_STRING(filename, "the name of the file with the tree topology data")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function constructs the hierarchical tree and all (interconnected) nodes forming the
        tree as described for the corresponding pure virtual function in the base class. For this
        class, this function recontructs the tree described in the configured tree topology data
        file. */
    vector<TreeNode*> constructTree() override;
};

//////////////////////////////////////////////////////////////////////

#endif
