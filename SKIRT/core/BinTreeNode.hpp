/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef BINTREENODE_HPP
#define BINTREENODE_HPP

#include "TreeNode.hpp"

//////////////////////////////////////////////////////////////////////

/** BinTreeNode is a TreeNode subclass that represents nodes in a binary tree, where each
    subsequent level subdivides a node into two equal portions along one of the coordinate planes.
    The class implements an alternating subdivision scheme, i.e. the nodes are alternatively
    divided perpendicular to each of the three coordinate axes when descending the tree. */
class BinTreeNode : public TreeNode
{
public:
    /** This constructor creates a new binary tree node with the specified parent node, identifier,
        and spatial extent. The constructor sets the level of the new node to be one higher than
        the level of the parent. If the pointer to the parent is null, the level of the new node is
        zero. */
    using TreeNode::TreeNode;

public:
    /** This function creates two new nodes subdividing the node at its geometric center along a
        plane perpendicular to one of the coordinate axes, depending on the node's level in the
        tree. The splitting direction is selected as the modulo of the node's level, with (0=x,
        1=y, 2=z), so that the nodes are alternatively divided along each of the axes when
        descending the tree.

        The children are assigned consecutive integer identifiers, starting with the identifier
        specified as an argument to this function. A node does NOT take ownership of its children,
        so the caller is responsible for deleting the child nodes when they are no longer needed.
        Invoking this function on a node that already has children results in undefined behavior.
        */
    void createChildren(int id) override;

    /** This function returns a pointer to the node's child that contains the specified point. More
        accurately, it returns the child corresponding to the half-space that contains the
        specified point relative to the node's central division plane. If the specified point is
        inside the node, then it will also be inside the returned child. Invoking this function on
        a childless node results in undefined behavior. */
    TreeNode* child(Vec r) override;

    /** This function adds the relevant neighbors to a node with children (the function does
        nothing if the node doesn't have any children). It considers internal neighbors among the
        children as well as the neighbors of the parent node (i.e. this node). These external
        neighbors are distributed among the children depending on the geometry; note that a
        particular external neighbor may be neighbor to multiple children. */
    void addNeighbors() override;
};

//////////////////////////////////////////////////////////////////////

#endif
