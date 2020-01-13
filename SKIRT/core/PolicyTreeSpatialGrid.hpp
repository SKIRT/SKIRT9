/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef POLICYTREESPATIALGRID_HPP
#define POLICYTREESPATIALGRID_HPP

#include "TreePolicy.hpp"
#include "TreeSpatialGrid.hpp"

//////////////////////////////////////////////////////////////////////

/** PolicyTreeSpatialGrid is a concrete subclass of the TreeSpatialGrid class representing
    three-dimensional spatial grids with cuboidal cells organized in a hierarchical tree. The
    tree's root node encloses the complete spatial domain, and nodes on subsequent levels
    recursively divide space into ever finer nodes. The depth of the tree can vary from place to
    place. The leaf nodes (those that are not further subdivided) are the actual spatial cells.

    This class offers a user-configurable option to select the tree type. Depending on the type of
    node employed, the tree can become an octtree (8 children per node) or a binary tree (2
    children per node). Other node types could be implemented, as long as they are cuboids lined up
    with the coordinate axes.

    The actual tree construction, other than the choice of node type, is delegated to a
    user-configurable \em policy instance. Encapsulating the configurable options and the
    corresponding implementation mechanisms for constructing spatial tree grids in a separate class
    hierarchy allows offering and possibly combining different policies without complicating the
    TreeSpatialGrid class hierarchy. */
class PolicyTreeSpatialGrid : public TreeSpatialGrid
{
    /** The enumeration type indicating the type of tree to be constructed: an octtree (8 children
        per node) or a binary tree (2 children per node). */
    ENUM_DEF(TreeType, OctTree, BinTree)
        ENUM_VAL(TreeType, OctTree, "an octtree (8 children per node)")
        ENUM_VAL(TreeType, BinTree, "a binary tree (2 children per node)")
    ENUM_END()

    ITEM_CONCRETE(PolicyTreeSpatialGrid, TreeSpatialGrid, "a tree-based spatial grid")

        PROPERTY_ENUM(treeType, TreeType, "the type of tree")
        ATTRIBUTE_DEFAULT_VALUE(treeType, "OctTree")
        ATTRIBUTE_DISPLAYED_IF(treeType, "Level2")
        ATTRIBUTE_INSERT(treeType, "treeTypeOctTree:OctTreeGrid;treeTypeBinTree:BinTreeGrid")

        PROPERTY_ITEM(policy, TreePolicy, "the tree construction policy (configuration options)")
        ATTRIBUTE_DEFAULT_VALUE(policy, "DensityTreePolicy")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function constructs the hierarchical tree and all (interconnected) nodes forming the
        tree as described for the corresponding pure virtual function in the base class. For this
        class, this function merely creates a root node of the appropriate type for the tree type
        configured by the user, and then invokes the constructTree() function of the \em policy
        configured by the user. */
    vector<TreeNode*> constructTree() override;
};

//////////////////////////////////////////////////////////////////////

#endif
