/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "DensityTreePolicy.hpp"
#include "Log.hpp"
#include "TreeNode.hpp"

////////////////////////////////////////////////////////////////////

namespace
{
    const int maxLevel = 8;

    bool needsSubdivide(TreeNode* node)
    {
        // TO DO
        return node->level() < maxLevel;
    }
}

////////////////////////////////////////////////////////////////////

vector<TreeNode*> DensityTreePolicy::constructTree(TreeNode* root)
{
    // initialize the tree node list with the root node as the first item
    vector<TreeNode*> nodev{root};

    // recursively subdivide the root node until all nodes satisfy the configured criteria
    Log* log = find<Log>();
    log->infoSetElapsed(0);
    int currentlevel = -1;
    for (size_t l = 0; l < nodev.size(); ++l)   // note that nodev.size() changes along the way
    {
        if (l%5000 == 0) log->infoIfElapsed("Working on node number " + std::to_string(l), 0);

        TreeNode* node = nodev[l];
        if (node->isChildless() && needsSubdivide(node))
        {
            int level = node->level();
            if (level>currentlevel)
            {
                log->info("Starting subdivision to level " + std::to_string(level+1));
                currentlevel = level;
            }

            node->createChildren(nodev.size());
            node->addNeighbors();
            nodev.insert(nodev.end(), node->children().begin(), node->children().end());
        }
    }

    // sort the neighbors for all nodes
    for (auto node : nodev) node->sortNeighbors();
    return nodev;
}

////////////////////////////////////////////////////////////////////
