/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "TreeSpatialGrid.hpp"
#include "Log.hpp"
#include "PathSegmentGenerator.hpp"
#include "Random.hpp"
#include "SpatialGridPath.hpp"
#include "SpatialGridPlotFile.hpp"
#include "StringUtils.hpp"
#include "TreeNode.hpp"

////////////////////////////////////////////////////////////////////

TreeSpatialGrid::~TreeSpatialGrid()
{
    for (auto node : _nodev) delete node;
}

////////////////////////////////////////////////////////////////////

void TreeSpatialGrid::setupSelfAfter()
{
    BoxSpatialGrid::setupSelfAfter();

    // determine a small fraction relative to the spatial extent of the grid; used during path traversal
    _eps = 1e-12 * extent().widths().norm();

    // make subclass construct the tree
    Log* log = find<Log>();
    log->info("Constructing the spatial tree grid...");
    _nodev = constructTree();

    // construct the vectors to help translating between node indices (leaf and nonleaf) and cell indices (leaf only)
    //  _cellindexv : cell index m corresponding to each node in nodev; -1 for nonleaf nodes
    //  _idv;       : index in nodev for each cell (i.e. leaf node); corresponds to node ID
    int numNodes = _nodev.size();
    int m = 0;
    _cellindexv.resize(numNodes, -1);
    for (int l = 0; l != numNodes; ++l)
    {
        if (_nodev[l]->isChildless())
        {
            _idv.push_back(l);
            _cellindexv[l] = m;
            m++;
        }
    }

    // determine the number of cells at each level in the tree hierarchy
    vector<int> countv;
    int numCells = _idv.size();
    for (int m = 0; m != numCells; ++m)
    {
        int level = nodeForCellIndex(m)->level();
        if (level + 1 > static_cast<int>(countv.size())) countv.resize(level + 1);
        countv[level]++;
    }

    // log these statistics, including a basic histogram
    log->info("Finished construction of the spatial tree grid");
    log->info("Number of cells at each level in the tree hierarchy:");
    int numLevels = countv.size();
    int maxCount = *std::max_element(countv.cbegin(), countv.cend());
    for (int level = 0; level != numLevels; ++level)
    {
        size_t numStars = std::round(20. * countv[level] / maxCount);
        log->info("  Level " + StringUtils::toString(level, 'd', 0, 2) + ":"
                  + StringUtils::toString(countv[level], 'd', 0, 9) + " ("
                  + StringUtils::toString(100. * countv[level] / numCells, 'f', 1, 5) + "%)  |"
                  + string(numStars, '*'));
    }
    log->info("  TOTAL   :" + StringUtils::toString(numCells, 'd', 0, 9) + " (100.0%)");
}

////////////////////////////////////////////////////////////////////

int TreeSpatialGrid::numCells() const
{
    return _idv.size();
}

////////////////////////////////////////////////////////////////////

double TreeSpatialGrid::volume(int m) const
{
    return nodeForCellIndex(m)->volume();
}

////////////////////////////////////////////////////////////////////

double TreeSpatialGrid::diagonal(int m) const
{
    return nodeForCellIndex(m)->diagonal();
}

////////////////////////////////////////////////////////////////////

int TreeSpatialGrid::cellIndex(Position bfr) const
{
    const TreeNode* node = root()->leafChild(bfr);
    return node ? cellIndexForNode(node) : -1;
}

////////////////////////////////////////////////////////////////////

Position TreeSpatialGrid::centralPositionInCell(int m) const
{
    return Position(nodeForCellIndex(m)->extent().center());
}

////////////////////////////////////////////////////////////////////

Position TreeSpatialGrid::randomPositionInCell(int m) const
{
    return random()->position(nodeForCellIndex(m)->extent());
}

//////////////////////////////////////////////////////////////////////

class TreeSpatialGrid::MySegmentGenerator : public PathSegmentGenerator
{
    const TreeSpatialGrid* _grid{nullptr};
    const TreeNode* _node{nullptr};

public:
    MySegmentGenerator(const TreeSpatialGrid* grid) : _grid(grid) {}

    bool next() override
    {
        switch (state())
        {
            case State::Unknown:
            {
                // try moving the photon packet inside the grid; if this is impossible, return an empty path
                if (!moveInside(_grid->extent(), _grid->_eps)) return false;

                // get the node containing the current location;
                _node = _grid->root()->leafChild(r());

                // if the photon packet started outside the grid, return the corresponding nonzero-length segment;
                // otherwise fall through to determine the first actual segment
                if (ds() > 0.) return true;
            }

            // intentionally falls through
            case State::Inside:
            {
                // determine the segment from the current position to the first cell wall
                // and adjust the position and cell indices accordingly
                double xnext = (kx() < 0.0) ? _node->xmin() : _node->xmax();
                double ynext = (ky() < 0.0) ? _node->ymin() : _node->ymax();
                double znext = (kz() < 0.0) ? _node->zmin() : _node->zmax();
                double dsx = (fabs(kx()) > 1e-15) ? (xnext - rx()) / kx() : DBL_MAX;
                double dsy = (fabs(ky()) > 1e-15) ? (ynext - ry()) / ky() : DBL_MAX;
                double dsz = (fabs(kz()) > 1e-15) ? (znext - rz()) / kz() : DBL_MAX;

                double ds;
                TreeNode::Wall wall;
                if (dsx <= dsy && dsx <= dsz)
                {
                    ds = dsx;
                    wall = (kx() < 0.0) ? TreeNode::BACK : TreeNode::FRONT;
                }
                else if (dsy <= dsx && dsy <= dsz)
                {
                    ds = dsy;
                    wall = (ky() < 0.0) ? TreeNode::LEFT : TreeNode::RIGHT;
                }
                else
                {
                    ds = dsz;
                    wall = (kz() < 0.0) ? TreeNode::BOTTOM : TreeNode::TOP;
                }
                propagater(ds + _grid->_eps);
                setSegment(_grid->cellIndexForNode(_node), ds);

                // attempt to find the new node among the neighbors of the current node;
                // this should not fail unless the new location is outside the grid,
                // however on rare occasions it fails due to rounding errors (e.g. in a corner),
                // thus we use top-down search as a fall-back
                const TreeNode* oldnode = _node;
                _node = _node->neighbor(wall, r());
                if (!_node) _node = _grid->root()->leafChild(r());

                // if we're stuck in the same node,
                // try to escape by advancing the position to the next representable coordinates
                if (_node == oldnode)
                {
                    // try to escape by advancing the position to the next representable coordinates
                    propagateToNextAfter();
                    _node = _grid->root()->leafChild(r());
                }

                // if we're outside the domain or still stuck in the same node, terminate the path
                if (!_node || _node == oldnode) setState(State::Outside);
                return true;
            }

            case State::Outside:
            {
            }
        }
        return false;
    }
};

////////////////////////////////////////////////////////////////////

std::unique_ptr<PathSegmentGenerator> TreeSpatialGrid::createPathSegmentGenerator() const
{
    return std::make_unique<MySegmentGenerator>(this);
}

////////////////////////////////////////////////////////////////////

namespace
{
    // this function writes a "0" for a leaf node or a "1" for a nonleaf node
    // followed by the recursive topological representation of its children
    void writeTopologyForNode(TreeNode* node, TextOutFile* outfile)
    {
        if (node->isChildless())
            outfile->writeLine("0");
        else
        {
            outfile->writeLine("1");
            for (auto child : node->children()) writeTopologyForNode(child, outfile);
        }
    }
}

////////////////////////////////////////////////////////////////////

void TreeSpatialGrid::writeTopology(TextOutFile* outfile) const
{
    outfile->writeLine("# Topology for tree spatial grid with " + std::to_string(numCells()) + " cells");
    outfile->writeLine(std::to_string(root()->children().size()));  // zero if the root node is not subdivided
    writeTopologyForNode(root(), outfile);
}

////////////////////////////////////////////////////////////////////

void TreeSpatialGrid::write_xy(SpatialGridPlotFile* outfile) const
{
    // Output the root cell and all leaf cells that are close to the section plane
    outfile->writeRectangle(xmin(), ymin(), xmax(), ymax());
    int nCells = numCells();
    for (int m = 0; m != nCells; ++m)
    {
        TreeNode* node = nodeForCellIndex(m);
        if (fabs(node->zmin()) < 1e-8 * extent().zwidth())
        {
            outfile->writeRectangle(node->xmin(), node->ymin(), node->xmax(), node->ymax());
        }
    }
}

////////////////////////////////////////////////////////////////////

void TreeSpatialGrid::write_xz(SpatialGridPlotFile* outfile) const
{
    // Output the root cell and all leaf cells that are close to the section plane
    outfile->writeRectangle(xmin(), zmin(), xmax(), zmax());
    int nCells = numCells();
    for (int m = 0; m != nCells; ++m)
    {
        TreeNode* node = nodeForCellIndex(m);
        if (fabs(node->ymin()) < 1e-8 * extent().ywidth())
        {
            outfile->writeRectangle(node->xmin(), node->zmin(), node->xmax(), node->zmax());
        }
    }
}

////////////////////////////////////////////////////////////////////

void TreeSpatialGrid::write_yz(SpatialGridPlotFile* outfile) const
{
    // Output the root cell and all leaf cells that are close to the section plane
    outfile->writeRectangle(ymin(), zmin(), ymax(), zmax());
    int nCells = numCells();
    for (int m = 0; m != nCells; ++m)
    {
        TreeNode* node = nodeForCellIndex(m);
        if (fabs(node->xmin()) < 1e-8 * extent().xwidth())
        {
            outfile->writeRectangle(node->ymin(), node->zmin(), node->ymax(), node->zmax());
        }
    }
}

////////////////////////////////////////////////////////////////////

void TreeSpatialGrid::write_xyz(SpatialGridPlotFile* outfile) const
{
    // determine the number of cells at each level in the tree hierarchy
    vector<int> countv;
    int nCells = numCells();
    for (int m = 0; m != nCells; ++m)
    {
        int level = nodeForCellIndex(m)->level();
        if (level + 1 > static_cast<int>(countv.size())) countv.resize(level + 1);
        countv[level]++;
    }

    // determine the number of levels to be included in output
    int maxLevel = static_cast<int>(countv.size()) - 1;
    int highestWriteLevel = 0;
    int cumulativeCells = 0;
    for (; highestWriteLevel <= maxLevel; ++highestWriteLevel)
    {
        cumulativeCells += countv[highestWriteLevel];
        if (cumulativeCells > 2500) break;  // empirical number
    }

    // inform the user if we are limiting output
    if (highestWriteLevel < maxLevel)
        find<Log>()->info("Limiting 3D grid plot output tree to level " + std::to_string(highestWriteLevel) + ", i.e. "
                          + std::to_string(cumulativeCells) + " cells.");

    // output all leaf cells up to a certain level
    for (int m = 0; m != nCells; ++m)
    {
        TreeNode* node = nodeForCellIndex(m);
        if (node->level() <= highestWriteLevel)
            outfile->writeCube(node->xmin(), node->ymin(), node->zmin(), node->xmax(), node->ymax(), node->zmax());
    }
}

////////////////////////////////////////////////////////////////////

TreeNode* TreeSpatialGrid::root() const
{
    return _nodev[0];
}

////////////////////////////////////////////////////////////////////

TreeNode* TreeSpatialGrid::nodeForCellIndex(int m) const
{
    return _nodev[_idv[m]];
}

////////////////////////////////////////////////////////////////////

int TreeSpatialGrid::cellIndexForNode(const TreeNode* node) const
{
    return _cellindexv[node->id()];
}

////////////////////////////////////////////////////////////////////
