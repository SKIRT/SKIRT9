/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "PolicyTreeSpatialGrid.hpp"
#include "BinTreeNode.hpp"
#include "OctTreeNode.hpp"

////////////////////////////////////////////////////////////////////

vector<TreeNode*> PolicyTreeSpatialGrid::constructTree()
{
    // create the root node using the requested type
    TreeNode* root = nullptr;
    switch (_treeType)
    {
        case TreeType::OctTree: root = new OctTreeNode(extent()); break;
        case TreeType::BinTree: root = new BinTreeNode(extent()); break;
    }

    // tell policy to construct the tree
    return _policy->constructTree(root);
}

////////////////////////////////////////////////////////////////////
