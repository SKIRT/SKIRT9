/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "BoxSearch.hpp"
#include "EntityCollection.hpp"
#include "NR.hpp"

////////////////////////////////////////////////////////////////////

namespace
{
    // This function determines the grid block separation points along a specified axis 1->x, 2->y, 3->z
    // to ensure that the entities are evenly distributed across the grid blocks.
    // It does this by binning the bounding box centers at a high resolution and then calculating
    // the cumulative distribution to set the grid separation points.
    void makegrid(Array& grid, const vector<Box>& boxv, int axis, int gridsize, double cmin, double cmax)
    {
        // guard against degenerate zero-width entities that are all lined up along this coordinate
        if (cmin == cmax)
        {
            double eps = 1e-12 * (cmin ? abs(cmin) : 1.);
            cmin -= eps;
            cmax += eps;
        }

        // determine the entity center distribution by binning at a decent resolution
        int nbins = gridsize * 100;
        double binwidth = (cmax - cmin) / nbins;
        vector<int> bins(nbins);
        for (const auto& box : boxv)
        {
            double center = 0.;
            switch (axis)
            {
                case 1: center = box.center().x(); break;
                case 2: center = box.center().y(); break;
                case 3: center = box.center().z(); break;
            }
            bins[static_cast<int>((center - cmin) / binwidth)] += 1;
        }

        // determine grid separation points based on the cumulative distribution
        grid.resize(gridsize + 1);
        grid[0] = -std::numeric_limits<double>::infinity();
        int perblock = boxv.size() / gridsize;  // target number of entities per block
        int cumul = 0;                          // cumulative number of entities in processed bins
        int gridindex = 1;                      // index of the next grid separation point to be filled
        for (int binindex = 0; binindex < nbins; binindex++)
        {
            cumul += bins[binindex];
            if (cumul > perblock * gridindex)
            {
                grid[gridindex] = cmin + (binindex + 1) * binwidth;
                gridindex += 1;
                if (gridindex >= gridsize) break;
            }
        }
        grid[gridsize] = std::numeric_limits<double>::infinity();
    }
}

////////////////////////////////////////////////////////////////////

int BoxSearch::blockIndex(int i, int j, int k) const
{
    return ((i * _numBlocks) + j) * _numBlocks + k;
}

////////////////////////////////////////////////////////////////////

int BoxSearch::blockIndex(Vec bfr) const
{
    // get the indices for the block containing the position in each direction
    int i = NR::locateClip(_xgrid, bfr.x());
    int j = NR::locateClip(_ygrid, bfr.y());
    int k = NR::locateClip(_zgrid, bfr.z());

    // return the corresponding linear index
    return blockIndex(i, j, k);
}

////////////////////////////////////////////////////////////////////

BoxSearch::BoxSearch() {}

////////////////////////////////////////////////////////////////////

void BoxSearch::loadEntities(int numEntities, std::function<Box(int)> bounds,
                             std::function<bool(int, const Box&)> intersects)
{
    // abort if there are no entities
    if (numEntities <= 0)
    {
        _extent = Box();
        _numBlocks = 0;
        _minEntitiesPerBlock = 0;
        _maxEntitiesPerBlock = 0;
        _avgEntitiesPerBlock = 0;
        return;
    }

    // cache the bounding boxes because we need them a few times
    vector<Box> boxv;
    boxv.reserve(numEntities);
    for (int m = 0; m != numEntities; ++m) boxv.emplace_back(bounds(m));

    // calculate the extent of the search domain
    _extent = boxv[0];
    for (const auto& box : boxv) _extent.extend(box);

    // determine the number of blocks in each spatial direction;
    // the floor of 20 is exceeded only for more than a million entities
    _numBlocks = max(20, static_cast<int>(std::cbrt(numEntities) / 5.));

    // build the grids in each spatial direction
    makegrid(_xgrid, boxv, 1, _numBlocks, _extent.xmin(), _extent.xmax());
    makegrid(_ygrid, boxv, 2, _numBlocks, _extent.ymin(), _extent.ymax());
    makegrid(_zgrid, boxv, 3, _numBlocks, _extent.zmin(), _extent.zmax());

    // construct an empty entity list for each of the nb*nb*nb blocks
    _listv.clear();  // remove any pre-existing lists
    _listv.resize(_numBlocks * _numBlocks * _numBlocks);

    // add each entity to the list for every block that its bounding box overlaps
    for (int m = 0; m != numEntities; ++m)
    {
        const auto& box = boxv[m];

        // find indices for first and last grid block overlapped by bounding box, in each spatial direction
        int i1 = NR::locateClip(_xgrid, box.xmin());
        int i2 = NR::locateClip(_xgrid, box.xmax());
        int j1 = NR::locateClip(_ygrid, box.ymin());
        int j2 = NR::locateClip(_ygrid, box.ymax());
        int k1 = NR::locateClip(_zgrid, box.zmin());
        int k2 = NR::locateClip(_zgrid, box.zmax());

        // loop over all blocks in that 3D range
        for (int i = i1; i <= i2; i++)
            for (int j = j1; j <= j2; j++)
                for (int k = k1; k <= k2; k++)
                {
                    // add the entity to the list if it indeed overlaps the block
                    Box block(_xgrid[i], _ygrid[j], _zgrid[k], _xgrid[i + 1], _ygrid[j + 1], _zgrid[k + 1]);
                    if (intersects(m, block)) _listv[blockIndex(i, j, k)].push_back(m);
                }
    }

    // calculate statistics
    _minEntitiesPerBlock = numEntities;
    _maxEntitiesPerBlock = 0;
    int totalEntityReferences = 0;
    for (const auto& list : _listv)
    {
        int size = list.size();
        _minEntitiesPerBlock = min(_minEntitiesPerBlock, size);
        _maxEntitiesPerBlock = max(_maxEntitiesPerBlock, size);
        totalEntityReferences += size;
    }
    _avgEntitiesPerBlock = static_cast<double>(totalEntityReferences) / _listv.size();
}

////////////////////////////////////////////////////////////////////

const Box& BoxSearch::extent() const
{
    return _extent;
}

////////////////////////////////////////////////////////////////////

int BoxSearch::numBlocks() const
{
    return _numBlocks;
}

////////////////////////////////////////////////////////////////////

int BoxSearch::minEntitiesPerBlock() const
{
    return _minEntitiesPerBlock;
}

////////////////////////////////////////////////////////////////////

int BoxSearch::maxEntitiesPerBlock() const
{
    return _maxEntitiesPerBlock;
}

////////////////////////////////////////////////////////////////////

double BoxSearch::avgEntitiesPerBlock() const
{
    return _avgEntitiesPerBlock;
}

////////////////////////////////////////////////////////////////////

int BoxSearch::firstEntity(Vec bfr, std::function<bool(int)> contains) const
{
    if (_numBlocks)
    {
        // loop over all entities overlapping the block containing the position
        for (int m : _listv[blockIndex(bfr)])
        {
            // return the first entity that actually contains the position
            if (contains(m)) return m;
        }
    }

    // there is no entity containing the position
    return -1;
}

////////////////////////////////////////////////////////////////////

double BoxSearch::accumulate(Vec bfr, std::function<double(int)> weight) const
{
    double sum = 0.;

    if (_numBlocks)
    {
        // sum the weights for all entities overlapping the block containing the position
        for (int m : _listv[blockIndex(bfr)]) sum += weight(m);
    }
    return sum;
}

////////////////////////////////////////////////////////////////////

void BoxSearch::getEntities(EntityCollection& entities, Vec bfr, std::function<double(int)> weight) const
{
    entities.clear();

    if (_numBlocks)
    {
        // add all entities overlapping that block to the collection with their respective weight
        for (int m : _listv[blockIndex(bfr)]) entities.add(m, weight(m));
    }
}

////////////////////////////////////////////////////////////////////

void BoxSearch::getEntities(EntityCollection& entities, Position bfr, Direction bfk,
                            std::function<double(int)> weight) const
{
    entities.clear();

    if (_numBlocks)
    {
        // verify that the path intersects the domain
        double smin, smax;
        if (_extent.intersects(bfr, bfk, smin, smax))
        {
            // find the indices for first and last block, in each spatial direction,
            // overlapped by the bounding box of the path's intersection with the domain
            Box pathbox(bfr + smin * bfk, bfr + smax * bfk);
            int i1 = NR::locateClip(_xgrid, pathbox.xmin());
            int i2 = NR::locateClip(_xgrid, pathbox.xmax());
            int j1 = NR::locateClip(_ygrid, pathbox.ymin());
            int j2 = NR::locateClip(_ygrid, pathbox.ymax());
            int k1 = NR::locateClip(_zgrid, pathbox.zmin());
            int k2 = NR::locateClip(_zgrid, pathbox.zmax());
            if (i1 > i2) std::swap(i1, i2);  // fix reversed coords for negative bfk components
            if (j1 > j2) std::swap(j1, j2);
            if (k1 > k2) std::swap(k1, k2);

            // loop over all blocks in that 3D range
            for (int i = i1; i <= i2; i++)
                for (int j = j1; j <= j2; j++)
                    for (int k = k1; k <= k2; k++)
                    {
                        // if the path intersects the block
                        Box block(_xgrid[i], _ygrid[j], _zgrid[k], _xgrid[i + 1], _ygrid[j + 1], _zgrid[k + 1]);
                        if (block.intersects(bfr, bfk, smin, smax))
                        {
                            // add all entities overlapping that block to the collection with their respective weight
                            for (int m : _listv[blockIndex(i, j, k)])
                            {
                                entities.add(m, weight(m));
                            }
                        }
                    }
        }
    }
}

////////////////////////////////////////////////////////////////////
