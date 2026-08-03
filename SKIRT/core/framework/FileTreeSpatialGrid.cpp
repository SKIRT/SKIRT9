/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "FileTreeSpatialGrid.hpp"
#include "BinTreeNode.hpp"
#include "FatalError.hpp"
#include "FilePaths.hpp"
#include "Log.hpp"
#include "OctTreeNode.hpp"
#include "System.hpp"
#include <fstream>

////////////////////////////////////////////////////////////////////

namespace
{
    // if the next number in the input file is a "1", this function subdivides the specified node and then
    // recursively calls itself for all of the newly created children of that node;
    // if the next number in the input file is a "0", this function does nothing
    void subdivideNodeIfNeeded(TreeNode* node, std::ifstream& infile, vector<TreeNode*>& nodev, Log* log)
    {
        int needsSubdivision = -1;
        infile >> needsSubdivision;
        switch (needsSubdivision)
        {
            case 0: break;
            case 1:
                log->infoIfElapsed("Subdiving node " + std::to_string(node->id()), 0);
                node->subdivide(nodev);
                for (auto child : node->children()) subdivideNodeIfNeeded(child, infile, nodev, log);
                break;
            default: throw FATALERROR("Topology input file has improper format and/or missing data");
        }
    }
}

////////////////////////////////////////////////////////////////////

vector<TreeNode*> FileTreeSpatialGrid::constructTree()
{
    // open the input file
    string filepath = find<FilePaths>()->input(_filename);
    std::ifstream infile = System::ifstream(filepath);
    if (!infile) throw FATALERROR("Could not open the spatial tree grid topology text file " + filepath);

    // log "reading file" message
    auto log = find<Log>();
    log->info("Reading tree topology from text file " + filepath + "...");

    // skip any header lines
    string line;
    while (infile.peek() == '#') getline(infile, line);

    // create the root node using the appropriate type
    TreeNode* root = nullptr;
    int numChildren = -1;
    infile >> numChildren;
    switch (numChildren)
    {
        case 8:
        case 0: root = new OctTreeNode(extent()); break;
        case 2: root = new BinTreeNode(extent()); break;
        default: throw FATALERROR("Topology input file specifies unsupported number of children in node subdivision");
    }

    // initialize the tree node list with the root node as the first item
    vector<TreeNode*> nodev{root};

    // recursively subdivide nodes according to the topology described in the file
    log->infoSetElapsed(0);
    subdivideNodeIfNeeded(root, infile, nodev, log);

    // log "done reading file" message
    log->info("Done reading tree topology");
    return nodev;
}

////////////////////////////////////////////////////////////////////
