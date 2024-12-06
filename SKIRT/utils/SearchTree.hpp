
/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef SEARCHTREE_HPP
#define SEARCHTREE_HPP

#include "Box.hpp"
#include "Vec.hpp"
#include <functional>

/** The SearchTree class represents a \f$k\f$-d tree \f$(k=3)\f$ used for efficient spatial searches. 
    It organizes a set of 3-dimensional points, \em Vec, to allow for efficiently querying the nearest
    neighbor given some input \em Vec. This functionality is implemented in the \em nearest(Vec). For 
    more information, see the Wikipedia page on \f$k\f$-d trees: https://en.wikipedia.org/wiki/K-d_tree
    
    Each node in the tree is represented by the nested \em SearchTree::Node class. */
class SearchTree

{
public:
    // forward declaration of the Node class
    class Node;

    //================= Construction - Destruction - Moving =================

    /** Recursively destroys the search tree by deleting the root. */
    ~SearchTree();

    /** Constructor that initializes an empty SearchTree. It is not ready for queries until the tree is built
        using one of the buildTree functions. */
    SearchTree();

    /** Builds a SearchTree with the all node positions. */
    void buildTree(const vector<Vec>& nodePositions);

    /** Builds a SearchTree with the specified node positions. This will reorder the indices
        in the process of building the tree, the resulting order is not specified. */
    void buildTree(const vector<Vec>& nodePositions, vector<int>& indices);

private:
    /** Recursively builds the k-d tree from the specified range of points. The depth determines
        the axis along which the space is partitioned at each level of the tree. */
    Node* buildTree(const vector<Vec>& nodePositions, vector<int>::iterator first, vector<int>::iterator last,
                    int depth) const;

    //================= Interrogation =================

public:
    /** Returns the index that represents the node nearest to the query point. It recursively descends the 
        tree to find the nearest neighbor of the \em _nodePositions. */
    int nearest(Vec& bfr, const vector<Vec>& nodePositions) const;

    Node* root() const;

    //================= Node Class =================

public:
    /** The Node class represents a node in the \f$k\f$-d tree \f$(k=3)\f$. Each node stores an index to a point in
        the dataset, the axis along which the space is partitioned at that node, and pointers to
        its parent and two child nodes. */
    class Node
    {
    private:
        //================= Node: Construction - Destruction =================

        /** Constructor that initializes a node with the specified site index, depth, and child
            pointers. The depth determines the axis along which the space is partitioned. */
        Node(int m, int depth, Node* left, Node* right);

        /** Destructor that deletes the children of the node. */
        ~Node();

        /** Sets the parent pointer for the node. This function is called from the parent's
            constructor. */
        void setParent(Node* up);

        //================= Node: Interrogation =================

        /** Returns the site index stored in the node. */
        int m() const;

        /** Returns a pointer to the parent node. */
        Node* up() const;

        /** Returns a pointer to the left child node. */
        Node* left() const;

        /** Returns a pointer to the right child node. */
        Node* right() const;

        /** Returns the appropriate child node for the specified query point. The child is
            determined based on the partitioning axis and the position of the query point. */
        Node* child(Vec bfr, const vector<Vec>& nodePositions) const;

        /** Returns the other child node than the one that would be appropriate for the specified
            query point. This function is used during the search process to explore both sides of
            the partitioning plane. */
        Node* otherChild(Vec bfr, const vector<Vec>& nodePositions) const;

        /** Returns the positions of the site represented by this node. */
        const Vec& position(const vector<Vec>& nodePositions) const;

        /** Returns the squared distance from the query point to the partitioning plane at this
            node. This distance is used to determine whether to explore the other side of the
            partitioning plane during the search process. */
        double squaredDistanceToSplitPlane(Vec bfr, const vector<Vec>& nodePositions) const;

        /** Returns the node in this subtree that represents the site nearest to the query point.
            The search process recursively descends the tree, unwinding the recursion to find the
            nearest neighbor. */
        Node* nearest(Vec bfr, const vector<Vec>& nodePositions);

        //================= Node: Data Mebers =================

        int _m;        // index in the _nodePositions vector
        int _axis;     // split axis for this node (0,1,2)
        Node* _up;     // pointer to the parent node
        Node* _left;   // pointer to the left child node
        Node* _right;  // pointer to the right child node

        friend class SearchTree;  // only allow SearchTree to access private members
    };

    //================= Data Members =================

private:
    Node* _root{nullptr};  // pointer to the root node
};

#endif