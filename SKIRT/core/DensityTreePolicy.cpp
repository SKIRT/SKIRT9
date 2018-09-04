/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "DensityTreePolicy.hpp"
#include "Array.hpp"
#include "Log.hpp"
#include "Parallel.hpp"
#include "ParallelFactory.hpp"
#include "ProcessManager.hpp"
#include "TreeNode.hpp"

////////////////////////////////////////////////////////////////////

namespace
{
    // maximum number of nodes evaluated between two invocations of infoIfElapsed()
    const size_t logEvalChunkSize = 10000;

    // maximum number of nodes subdivided between two invocations of infoIfElapsed()
    const size_t logDivideChunkSize = 5000;
}

////////////////////////////////////////////////////////////////////

vector<TreeNode*> DensityTreePolicy::constructTree(TreeNode* root)
{
    auto log = find<Log>();
    auto parallel = find<ParallelFactory>()->parallelDistributed();

    // initialize the tree node list with the root node as the first item
    vector<TreeNode*> nodev{root};

    // recursively subdivide the root node until the minimum level has been reached
    int level = 0;      // current level
    size_t lbeg = 0;    // node index range for the current level;
    size_t lend = 1;    // at level 0, the node list contains just the root node
    while (level!=minLevel())
    {
        log->info("Subdividing level " + std::to_string(level) + ": " + std::to_string(lend-lbeg) + " nodes");
        log->infoSetElapsed(lend-lbeg);
        for (size_t l=lbeg; l!=lend; ++l)
        {
            nodev[l]->subdivide(nodev);
            if ((l+1)%logDivideChunkSize == 0)
                log->infoIfElapsed("Subdividing level " + std::to_string(level) + ": ", logDivideChunkSize);
        }
        // update iteration variables to the next level
        level++;
        lbeg = lend;
        lend = nodev.size();
    }

    // recursively subdivide the nodes beyond the minimum level until all nodes satisfy the configured criteria
    while (level!=maxLevel())
    {
        size_t numEvalNodes = lend-lbeg;
        log->info("Subdividing level " + std::to_string(level) + ": " + std::to_string(numEvalNodes) + " nodes");
        log->infoSetElapsed(numEvalNodes);

        // evaluate nodes at this level: value in the array becomes one for nodes that need to be subdivided
        // we parallelize this operation because it might be resource intensive (e.g. sampling densities)
        Array divide(numEvalNodes);
        parallel->call(numEvalNodes, [this, log, level, lbeg, &nodev, &divide](size_t firstIndex, size_t numIndices)
        {
            while (numIndices)
            {
                size_t currentChunkSize = min(logEvalChunkSize, numIndices);
                for (size_t l=firstIndex; l!=firstIndex+currentChunkSize; ++l)
                {
                    if (needsSubdivide(nodev[lbeg+l])) divide[l] = 1.;
                }
                log->infoIfElapsed("Evaluation for level " + std::to_string(level) + ": ", currentChunkSize);
                firstIndex += currentChunkSize;
                numIndices -= currentChunkSize;
            }
        });
        ProcessManager::sumToAll(divide);

        // subdivide the nodes that have been flagged
        size_t numDivideNodes = divide.sum();
        log->infoSetElapsed(numDivideNodes);
        size_t numDone = 0;
        for (size_t l=0; l!=numEvalNodes; ++l)
        {
            if (divide[l])
            {
                nodev[lbeg+l]->subdivide(nodev);
                numDone++;
                if (numDone%logDivideChunkSize == 0)
                    log->infoIfElapsed("Subdivision for level " + std::to_string(level) + ": ", logDivideChunkSize);
            }
        }
        // update iteration variables to the next level
        level++;
        lbeg = lend;
        lend = nodev.size();
    }

    // sort the neighbors for all nodes
    for (auto node : nodev) node->sortNeighbors();
    return nodev;
}

////////////////////////////////////////////////////////////////////

bool DensityTreePolicy::needsSubdivide(TreeNode* node)
{
    return (node->id() % 2) == 0;
}

////////////////////////////////////////////////////////////////////
