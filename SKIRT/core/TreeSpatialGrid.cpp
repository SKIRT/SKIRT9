/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "TreeSpatialGrid.hpp"
#include "Log.hpp"
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
    for (int l=0; l!=numNodes; ++l)
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
    for (int m=0; m!=numCells; ++m)
    {
        int level = nodeForCellIndex(m)->level();
        if (level+1 > static_cast<int>(countv.size())) countv.resize(level+1);
        countv[level]++;
    }

    // log these statistics
    log->info("Finished construction of the spatial tree grid");
    log->info("Number of cells at each level in the tree hierarchy:");
    int numLevels = countv.size();
    for (int level=0; level!=numLevels; ++level)
        log->info("  Level " + StringUtils::toString(level, 'd', 0, 2) + ":"
                             + StringUtils::toString(countv[level], 'd', 0, 9) + " ("
                             + StringUtils::toString(100.*countv[level]/numCells, 'f', 1, 5) + "%)");
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

////////////////////////////////////////////////////////////////////

void TreeSpatialGrid::path(SpatialGridPath* path) const
{
    // initialize the path
    path->clear();

    // if the photon package starts outside the dust grid, move it into the first grid cell that it will pass
    Position bfr = path->moveInside(extent(), _eps);

    // get the node containing the current location;
    // if the position is not inside the grid, return an empty path
    const TreeNode* node = root()->leafChild(bfr);
    if (!node) return path->clear();

    // get the starting point and direction
    double x,y,z;
    bfr.cartesian(x,y,z);
    double kx,ky,kz;
    path->direction().cartesian(kx,ky,kz);

    // loop over nodes/path segments until we leave the grid
    while (node)
    {
        double xnext = (kx<0.0) ? node->xmin() : node->xmax();
        double ynext = (ky<0.0) ? node->ymin() : node->ymax();
        double znext = (kz<0.0) ? node->zmin() : node->zmax();
        double dsx = (fabs(kx)>1e-15) ? (xnext-x)/kx : DBL_MAX;
        double dsy = (fabs(ky)>1e-15) ? (ynext-y)/ky : DBL_MAX;
        double dsz = (fabs(kz)>1e-15) ? (znext-z)/kz : DBL_MAX;

        double ds;
        TreeNode::Wall wall;
        if (dsx<=dsy && dsx<=dsz)
        {
            ds = dsx;
            wall = (kx<0.0) ? TreeNode::BACK : TreeNode::FRONT;
        }
        else if (dsy<=dsx && dsy<=dsz)
        {
            ds = dsy;
            wall = (ky<0.0) ? TreeNode::LEFT : TreeNode::RIGHT;
        }
        else
        {
            ds = dsz;
            wall = (kz<0.0) ? TreeNode::BOTTOM : TreeNode::TOP;
        }
        path->addSegment(cellIndexForNode(node), ds);
        x += (ds+_eps)*kx;
        y += (ds+_eps)*ky;
        z += (ds+_eps)*kz;

        // attempt to find the new node among the neighbors of the current node;
        // this should not fail unless the new location is outside the grid,
        // however on rare occasions it fails due to rounding errors (e.g. in a corner),
        // thus we use top-down search as a fall-back
        const TreeNode* oldnode = node;
        node = node->neighbor(wall, Vec(x,y,z));
        if (!node) node = root()->leafChild(Vec(x,y,z));

        // if we're stuck in the same node...
        if (node==oldnode)
        {
            // try to escape by advancing the position to the next representable coordinates
            find<Log>()->warning("Photon package seems stuck in spatial cell "
                                 + std::to_string(node->id()) + " -- escaping");
            x = std::nextafter(x, (kx<0.0) ? -DBL_MAX : DBL_MAX);
            y = std::nextafter(y, (ky<0.0) ? -DBL_MAX : DBL_MAX);
            z = std::nextafter(z, (kz<0.0) ? -DBL_MAX : DBL_MAX);
            node = root()->leafChild(Vec(x,y,z));

            // if that didn't work, terminate the path
            if (node==oldnode)
            {
                find<Log>()->warning("Photon package is stuck in spatial cell "
                                     + std::to_string(node->id()) + " -- terminating this path");
                break;
            }
        }
    }
}

////////////////////////////////////////////////////////////////////

void TreeSpatialGrid::write_xy(SpatialGridPlotFile* outfile) const
{
    // Output the root cell and all leaf cells that are close to the section plane
    outfile->writeRectangle(xmin(), ymin(), xmax(), ymax());
    int nCells = numCells();
    for (int m=0; m!=nCells; ++m)
    {
        TreeNode* node = nodeForCellIndex(m);
        if (fabs(node->zmin()) < 1e-8*extent().zwidth())
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
    for (int m=0; m!=nCells; ++m)
    {
        TreeNode* node = nodeForCellIndex(m);
        if (fabs(node->ymin()) < 1e-8*extent().ywidth())
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
    for (int m=0; m!=nCells; ++m)
    {
        TreeNode* node = nodeForCellIndex(m);
        if (fabs(node->xmin()) < 1e-8*extent().xwidth())
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
    for (int m=0; m!=nCells; ++m)
    {
        int level = nodeForCellIndex(m)->level();
        if (level+1 > static_cast<int>(countv.size())) countv.resize(level+1);
        countv[level]++;
    }

    // determine the number of levels to be included in output
    int maxLevel = static_cast<int>(countv.size()) - 1;
    int highestWriteLevel = 0;
    int cumulativeCells = 0;
    for (; highestWriteLevel<=maxLevel; ++highestWriteLevel)
    {
        cumulativeCells += countv[highestWriteLevel];
        if (cumulativeCells > 2500) break;          // empirical number
    }

    // inform the user if we are limiting output
    if (highestWriteLevel < maxLevel)
        find<Log>()->info("Limiting 3D grid plot output tree to level " + std::to_string(highestWriteLevel) +
                  ", i.e. " + std::to_string(cumulativeCells) + " cells.");

    // output all leaf cells up to a certain level
    for (int m=0; m!=nCells; ++m)
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
