/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef TREEPOLICY_HPP
#define TREEPOLICY_HPP

#include "SimulationItem.hpp"
class TreeNode;

//////////////////////////////////////////////////////////////////////

/** TreePolicy is an abstract class that represents the configurable options and the corresponding
    implementation mechanisms for constructing spatial tree grids, as a service to the
    PolicyTreeSpatialGrid class. Encapsulating these policies in a separate class hierarchy allows
    offering and possibly combining different policies without complicating the TreeSpatialGrid
    class hierarchy.

    The public interface of TreePolicy class offers just a single function that is invoked by the
    PolicyTreeSpatialGrid class during setup to construct the tree hierarchy and all TreeNode
    instances in it. */
class TreePolicy : public SimulationItem
{
    ITEM_ABSTRACT(TreePolicy, SimulationItem, "a spatial tree grid construction policy")

        PROPERTY_INT(minLevel, "the minimum level of grid refinement")
        ATTRIBUTE_MIN_VALUE(minLevel, "0")
        ATTRIBUTE_MAX_VALUE(minLevel, "99")
        ATTRIBUTE_DEFAULT_VALUE(minLevel, "OctTreeGrid:3;BinTreeGrid:9")

        PROPERTY_INT(maxLevel, "the maximum level of grid refinement")
        ATTRIBUTE_MIN_VALUE(maxLevel, "0")
        ATTRIBUTE_MAX_VALUE(maxLevel, "99")
        ATTRIBUTE_DEFAULT_VALUE(maxLevel, "OctTreeGrid:7;BinTreeGrid:21")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function verifies that the maximum level is not below the minimum level. */
    void setupSelfBefore() override;

public:
    /** This function must be implemented in a subclass. It constructs the hierarchical tree and
        all (interconnected) nodes forming the tree. All nodes are instances of the same TreeNode
        subclass as the root node passed as an argument. The function returns a list of pointers to
        all created nodes (leaf and nonleaf). Ownership of the nodes resides in this returned list
        (not in the node hierarchy itself) and is thus handed to the caller.

        Each (leaf and nonleaf) node is given an identifier (ID) corresponding to its index in the
        list. The first node in the list is the root node of tree, i.e. the node passed as an
        argument. By definition, the root node has an ID of zero.

        This function also causes the nodes to construct neighbor lists, interconnecting the nodes
        according to their spatial relationship. The neighbor lists are sorted so that the
        neighbors with the largest overlapping border area are listed first, increasing (on
        average) the probability of locating the correct neighbor early in the list. */
    virtual vector<TreeNode*> constructTree(TreeNode* root) = 0;
};

//////////////////////////////////////////////////////////////////////

#endif
