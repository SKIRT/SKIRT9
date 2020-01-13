/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef TREENODE_HPP
#define TREENODE_HPP

#include "Box.hpp"
#include <array>

////////////////////////////////////////////////////////////////////

/** TreeNode is an abstract class that represents nodes in a TreeSpatialGrid. It holds a node
    identifier, the spatial extent of the node (a cuboid lined up with the coordinate axes), and
    links to the parent, children and neighbors of the node. */
class TreeNode : public Box
{
    //============= Constructing and destructing =============

public:
    /** This constructor creates a new tree node with the specified parent node, identifier, and
        spatial extent. The constructor sets the level of the new node to be one higher than the
        level of the parent. If the pointer to the parent is null, the level of the new node is
        zero. */
    TreeNode(TreeNode* parent, int id, const Box& extent);

    /** This constructor creates a new root node with the specified spatial extent. The constructor
        sets both the node identifier and the level of the new root node to zero. */
    TreeNode(const Box& extent);

    /** Trivial virtual destructor. */
    virtual ~TreeNode();

    //============= Basic getters =============

public:
    /** This function returns a pointer to the parent of the node. */
    TreeNode* parent();

    /** This function returns the ID number of the node. */
    int id() const;

    /** This function returns the level of the node. */
    int level() const;

    //============= Getting children =============

public:
    /** This function returns true if the node has no children, or false if it has children. */
    bool isChildless() const;

    /** This function returns a list of pointers to the node's children. The list is either empty or
        of the size appropriate for the subclass (e.g. 2 for binary tree, 8 for octtree). */
    const vector<TreeNode*>& children() const;

protected:
    /** This function returns a pointer to the node's child with index \f$l\f$. The appropriate
        range for the index depends on the subclass (e.g. 0-1 for binary tree, 0-7 for octtree).
        Invoking this function on a childless node of with a child index out of the appropriate
        range results in undefined behavior. */
    TreeNode* childAt(int l);

public:
    /** This function returns a pointer to the node's child that contains the specified point,
        assuming that the point is inside the node. Invoking this function on a childless node
        results in undefined behavior. */
    virtual TreeNode* child(Vec r) = 0;

    /** This function returns a pointer to the deepest node in the descendent hierarchy of this
        node that contains the specified position, or the null pointer if the position is outside
        the node. It uses the child(r) function resursively to locate the appropriate node. */
    TreeNode* leafChild(Vec r);

    //============= Managing children =============

public:
    /** This function subdivides the node by creating the appropriate number of child subnodes
        (through the createChildren() function) and appending pointers to these children to the
        specified node list. The node identifiers of the child nodes are set so that the node
        identifier matches the index of the node in the node list. Finally, the neighbor lists are
        updated (through the addNeighbors() function). */
    void subdivide(vector<TreeNode*>& nodev);

    /** This function creates new nodes partitioning the node, and adds these new nodes as its own
        child nodes. Subdivision happens according to a fixed scheme determined by each subclass.
        The children are assigned consecutive integer identifiers, starting with the identifier
        specified as an argument to this function. A node does NOT take ownership of its children,
        so the caller is responsible for deleting the child nodes when they are no longer needed.
        Invoking this function on a node that already has children results in undefined behavior.
        */
    virtual void createChildren(int id) = 0;

protected:
    /** This function adds the specified child to the end of the child list. */
    void addChild(TreeNode* child);

    //============= Getting neighbors =============

public:
    /** This enum contains a constant for each of the six walls in a node. The x-coordinate
        increases from BACK to FRONT, the y-coordinate increases from LEFT to RIGHT, and the
        z-coordinate increases from BOTTOM to TOP. */
    enum Wall { BACK = 0, FRONT, LEFT, RIGHT, BOTTOM, TOP };

    /** This function returns a list of pointers to the node's neighbors at the given wall. The
        list may be empty, for example because neighbors are still being added. */
    const vector<TreeNode*>& neighbors(Wall wall) const;

    /** This function returns a pointer to the node just beyond a given wall that contains the
        specified position, or null if such a node can't be found by searching the neighbors of
        that wall. The function expects that the neighbors of the node have been added. */
    const TreeNode* neighbor(Wall wall, Vec r) const;

    //============= Managing neighbors =============

public:
    /** This function adds the relevant neighbors to a node with children (the function does
        nothing if the node doesn't have any children). It considers internal neighbors among the
        children as well as the neighbors of the parent node (i.e. this node). These external
        neighbors are distributed among the children depending on the geometry; note that a
        particular external neighbor may be neighbor to multiple children. */
    virtual void addNeighbors() = 0;

    /** This function adds a node to the list of neighbors corresponding to a given wall. */
    void addNeighbor(Wall wall, TreeNode* node);

    /** This function deletes a node from the list of neighbors corresponding to a given wall. */
    void deleteNeighbor(Wall wall, TreeNode* node);

    /** This static function makes the two specified nodes neighbors by adding node2 as a neighbor
        to node1 at wall1, and adding node1 as a neighbor to node2 at the complementing wall
        (back/front, left/right, bottom/top). */
    static void makeNeighbors(Wall wall1, TreeNode* node1, TreeNode* node2);

    /** This function sorts the neighbor lists for each wall of this node so that neighbors with a
        larger overlap area are listed first. The function should be called only after neighbors
        have been added for all nodes in the tree. */
    void sortNeighbors();

    //============= Data members =============

private:
    int _id{0};
    int _level{0};
    TreeNode* _parent{nullptr};
    vector<TreeNode*> _children;
    std::array<vector<TreeNode*>, 6> _neighbors;  // one empty neighbor list for each wall
};

////////////////////////////////////////////////////////////////////

#endif
