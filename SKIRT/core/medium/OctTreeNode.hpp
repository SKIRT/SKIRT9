/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef OCTTREENODE_HPP
#define OCTTREENODE_HPP

#include "TreeNode.hpp"

//////////////////////////////////////////////////////////////////////

/** OctTreeNode is a TreeNode subclass that represents nodes in an octtree, where each subsequent
    level subdivides a node into eight equal portions (i.e. using a regular geometric subdivision
    scheme). */
class OctTreeNode : public TreeNode
{
public:
    /** This constructor creates a new octtree node with the specified parent node, identifier, and
        spatial extent. The constructor sets the level of the new node to be one higher than the
        level of the parent. If the pointer to the parent is null, the level of the new node is
        zero. */
    using TreeNode::TreeNode;

    /** This function creates eight new equally sized subnodes subdividing the node at the
        geometric center, and adds these new nodes as its own child nodes. The children are
        assigned consecutive integer identifiers, starting with the identifier specified as an
        argument to this function. A node does NOT take ownership of its children, so the caller is
        responsible for deleting the child nodes when they are no longer needed. Invoking this
        function on a node that already has children results in undefined behavior. */
    void createChildren(int id) override;

    /** This function returns a pointer to the node's child that contains the specified point. More
        accurately, it returns the child corresponding to the octant that contains the specified
        point relative to the node's central division point. If the specified point is inside the
        node, then it will also be inside the returned child. This function crashes if the node is
        childless. */
    TreeNode* child(Vec r) override;

    /** This function adds the relevant neighbors to a node with children (the function does
        nothing if the node doesn't have any children). It considers internal neighbors (each of
        the 8 children has 3 neighbors among its siblings) as well as the neighbors of the parent
        node (i.e. this node). These external neighbors are distributed among the children
        depending on the geometry; note that a particular external neighbor may be neighbor to
        multiple children. */
    void addNeighbors() override;
};

//////////////////////////////////////////////////////////////////////

#endif
