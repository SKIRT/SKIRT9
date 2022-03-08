/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "SiteListTreePolicy.hpp"
#include "Log.hpp"
#include "MediumSystem.hpp"
#include "SiteListInterface.hpp"
#include "TreeNode.hpp"

////////////////////////////////////////////////////////////////////

namespace
{
    // maximum number of nodes subdivided between two invocations of infoIfElapsed()
    const size_t logDivideChunkSize = 5000;

    // maximum number of sites inserted between two invocations of infoIfElapsed()
    const size_t logInsertChunkSize = 10000;
}

////////////////////////////////////////////////////////////////////

namespace
{
    // private function to insert a site into the tree, performing recursive subdivision if needed
    void insertSite(const SiteListInterface* sli, int newSite, TreeNode* parent, int maxLevel, vector<int>& sitev,
                    vector<TreeNode*>& nodev)
    {
        // find the leaf node that contains this site's position
        TreeNode* node = parent->leafChild(sli->sitePosition(newSite));
        if (node)
        {
            int id = node->id();

            // if the leaf node is still empty, just add the site to it
            if (sitev[id] < 0)
            {
                sitev[id] = newSite;
            }

            // if the leaf node already holds a site, and it is not at the maximum level,
            // subdivide the node, and insert both the old and new particles into the node
            else if (node->level() < maxLevel)
            {
                node->subdivide(nodev);
                int numChildNodes = node->children().size();
                while (numChildNodes--) sitev.push_back(-1);
                insertSite(sli, sitev[id], node, maxLevel, sitev, nodev);
                insertSite(sli, newSite, node, maxLevel, sitev, nodev);
            }
        }
    }
}

////////////////////////////////////////////////////////////////////

vector<TreeNode*> SiteListTreePolicy::constructTree(TreeNode* root)
{
    // locate the medium offering the site list
    auto sli = find<MediumSystem>()->interface<SiteListInterface>(2);
    auto log = find<Log>();

    // initialize the tree node list with the root node as the first item
    vector<TreeNode*> nodev{root};

    // recursively subdivide the root node until the minimum level has been reached
    size_t lbeg = 0;  // node index range for the current level;
    size_t lend = 1;  // at level 0, the node list contains just the root node
    for (int level = 0; level != minLevel(); ++level)
    {
        log->info("Subdividing level " + std::to_string(level) + ": " + std::to_string(lend - lbeg) + " nodes");
        log->infoSetElapsed(lend - lbeg);
        for (size_t l = lbeg; l != lend; ++l)
        {
            nodev[l]->subdivide(nodev);
            if ((l + 1) % logDivideChunkSize == 0)
                log->infoIfElapsed("Subdividing level " + std::to_string(level) + ": ", logDivideChunkSize);
        }
        // update iteration variables to the next level
        lbeg = lend;
        lend = nodev.size();
    }

    // initialize a list, used only during construction, that contains an integer value corresponding to
    // each node (nonleaf and leaf) so far created; for a nonleaf node the value is undefined;
    // for a leaf node that "holds" one of the sites, the value is the index of the site in the site list;
    // for a leaf node that does not yet "hold" a site, the value is -1
    vector<int> sitev(nodev.size(), -1);

    // add sites one by one, subdividing recursively if the leaf node containing the new site position
    // already "holds" another particle
    int numSites = sli->numSites();
    log->info("Subdividing tree to insert " + std::to_string(numSites) + " sites");
    log->infoSetElapsed(numSites);
    for (int i = 0; i < numSites; i++)
    {
        insertSite(sli, i, root, maxLevel(), sitev, nodev);
        if ((i + 1) % logInsertChunkSize == 0)
            log->infoIfElapsed("Inserting site " + std::to_string(i) + ": ", logInsertChunkSize);
    }

    // perform additional subdivisions as requested
    lbeg = 0;
    lend = nodev.size();
    for (int extra = 0; extra != numExtraLevels(); ++extra)
    {
        log->info("Subdividing extra level " + std::to_string(extra) + ": " + std::to_string(lend - lbeg) + " nodes");
        log->infoSetElapsed(lend - lbeg);
        for (size_t l = lbeg; l != lend; ++l)
        {
            if (nodev[l]->isChildless() && nodev[l]->level() < maxLevel()) nodev[l]->subdivide(nodev);
            if ((l + 1) % logDivideChunkSize == 0)
                log->infoIfElapsed("Subdividing extra level " + std::to_string(extra) + ": ", logDivideChunkSize);
        }
        // update iteration variables to the next level
        lbeg = lend;
        lend = nodev.size();
    }

    // sort the neighbors for all nodes
    for (auto node : nodev) node->sortNeighbors();
    return nodev;
}

////////////////////////////////////////////////////////////////////
