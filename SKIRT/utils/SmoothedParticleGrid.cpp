/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "SmoothedParticleGrid.hpp"
#include "NR.hpp"
#include "SmoothedParticle.hpp"
#include <set>

////////////////////////////////////////////////////////////////////

namespace
{
    // returns the linear index for cell (i,j,k) in a m*m*m table
    inline int index(int m, int i, int j, int k) { return ((i * m) + j) * m + k; }

    // builds a smart grid in the specified spatial direction (1=x, 2=y, 3=z) and with the specified size,
    // and stores it in output parameter "grid", along with the minimum and maximum coordinate enclosing the particles
    void makegrid(const vector<SmoothedParticle>& pv, int dir, int gridsize, Array& grid, double& cmin, double& cmax)
    {
        int n = pv.size();

        // find the spatial range of the particles in the specified direction
        cmin = +std::numeric_limits<double>::infinity();
        cmax = -std::numeric_limits<double>::infinity();
        for (int p = 0; p < n; p++)
        {
            cmin = min(cmin, pv[p].center(dir) - pv[p].radius());
            cmax = max(cmax, pv[p].center(dir) + pv[p].radius());
        }

        // guard against point sources (h=0) that are all lined up along this coordinate
        if (cmin == cmax)
        {
            double eps = 1e-12 * (cmin ? abs(cmin) : 1.);
            cmin -= eps;
            cmax += eps;
        }

        // determine the particle distribution by binning at a decent resolution
        int nbins = gridsize * 100;
        double binwidth = (cmax - cmin) / nbins;
        vector<int> bins(nbins);
        for (int p = 0; p < n; p++) bins[int((pv[p].center(dir) - cmin) / binwidth)] += 1;

        // determine grid separation points based on the cumulative distribution
        grid.resize(gridsize + 1);
        grid[0] = -std::numeric_limits<double>::infinity();
        int percell = n / gridsize;  // target number of particles per cell
        int cumul = 0;               // cumulative number of particles in processed bins
        int gridindex = 1;           // index of the next grid separation point to be filled
        for (int binindex = 0; binindex < nbins; binindex++)
        {
            cumul += bins[binindex];
            if (cumul > percell * gridindex)
            {
                grid[gridindex] = cmin + (binindex + 1) * binwidth;
                gridindex += 1;
                if (gridindex >= gridsize) break;
            }
        }
        grid[gridsize] = std::numeric_limits<double>::infinity();
    }

    // returns the square of the argument
    inline double square(double value) { return value * value; }

    // Determines whether an axis-aligned bounding box intersects with a sphere
    // (algorithm due to Jim Arvo in "Graphics Gems" (1990))
    // xmin, xmax, ymin, ymax, zmin, zmax: ll and ur corner of the bounding box
    // xc, yc, zc, r: center and radius of the sphere
    // returns true if the bounding box and the sphere intersect, false otherwise
    bool intersects(double xmin, double xmax, double ymin, double ymax, double zmin, double zmax, double xc, double yc,
                    double zc, double r)
    {
        double squaredist = square(r);

        if (xc < xmin)
            squaredist -= square(xc - xmin);
        else if (xc > xmax)
            squaredist -= square(xc - xmax);
        if (yc < ymin)
            squaredist -= square(yc - ymin);
        else if (yc > ymax)
            squaredist -= square(yc - ymax);
        if (zc < zmin)
            squaredist -= square(zc - zmin);
        else if (zc > zmax)
            squaredist -= square(zc - zmax);

        return squaredist >= 0.;
    }
}

////////////////////////////////////////////////////////////////////

SmoothedParticleGrid::SmoothedParticleGrid(const vector<SmoothedParticle>& pv, int gridsize) : _m(gridsize)
{
    // build the grids in each spatial direction
    double xmin, ymin, zmin, xmax, ymax, zmax;
    makegrid(pv, 1, gridsize, _xgrid, xmin, xmax);
    makegrid(pv, 2, gridsize, _ygrid, ymin, ymax);
    makegrid(pv, 3, gridsize, _zgrid, zmin, zmax);
    setExtent(xmin, ymin, zmin, xmax, ymax, zmax);

    // make room for m*m*m cells
    _listv.resize(gridsize * gridsize * gridsize);

    // add each particle to the list for every cell that it overlaps
    int n = pv.size();
    for (int p = 0; p < n; p++)
    {
        Vec rc = pv[p].center();
        double h = pv[p].radius();

        // find indices for first and last cell possibly overlapped by particle, in each spatial direction
        int i1 = NR::locateClip(_xgrid, rc.x() - h);
        int i2 = NR::locateClip(_xgrid, rc.x() + h);
        int j1 = NR::locateClip(_ygrid, rc.y() - h);
        int j2 = NR::locateClip(_ygrid, rc.y() + h);
        int k1 = NR::locateClip(_zgrid, rc.z() - h);
        int k2 = NR::locateClip(_zgrid, rc.z() + h);

        // loop over all cells in that 3D range
        for (int i = i1; i <= i2; i++)
            for (int j = j1; j <= j2; j++)
                for (int k = k1; k <= k2; k++)
                {
                    // add the particle to the list if it indeed overlaps the cell
                    if (intersects(_xgrid[i], _xgrid[i + 1], _ygrid[j], _ygrid[j + 1], _zgrid[k], _zgrid[k + 1], rc.x(),
                                   rc.y(), rc.z(), h))
                        _listv[index(gridsize, i, j, k)].push_back(&pv[p]);
                }
    }

    // calculate statistics
    _pmin = n;
    _pmax = 0;
    _ptotal = 0;
    for (int index = 0; index < gridsize * gridsize * gridsize; index++)
    {
        int size = _listv[index].size();
        _pmin = min(_pmin, size);
        _pmax = max(_pmax, size);
        _ptotal += size;
    }
}

////////////////////////////////////////////////////////////////////

int SmoothedParticleGrid::minParticlesPerCell() const
{
    return _pmin;
}

////////////////////////////////////////////////////////////////////

int SmoothedParticleGrid::maxParticlesPerCell() const
{
    return _pmax;
}

////////////////////////////////////////////////////////////////////

int SmoothedParticleGrid::totalParticles() const
{
    return _ptotal;
}

////////////////////////////////////////////////////////////////////

const vector<const SmoothedParticle*>& SmoothedParticleGrid::particlesFor(Vec r) const
{
    int i = NR::locateClip(_xgrid, r.x());
    int j = NR::locateClip(_ygrid, r.y());
    int k = NR::locateClip(_zgrid, r.z());
    return _listv[index(_m, i, j, k)];
}

////////////////////////////////////////////////////////////////////

vector<const SmoothedParticle*> SmoothedParticleGrid::particlesFor(const Box& box) const
{
    // find indices for first and last cell possibly overlapping the box
    int i1 = NR::locateClip(_xgrid, box.xmin());
    int j1 = NR::locateClip(_ygrid, box.ymin());
    int k1 = NR::locateClip(_zgrid, box.zmin());
    int i2 = NR::locateClip(_xgrid, box.xmax());
    int j2 = NR::locateClip(_ygrid, box.ymax());
    int k2 = NR::locateClip(_zgrid, box.zmax());

    // if the box is fully inside a single cell, just return the corresponding list
    if (i1 == i2 && j1 == j2 && k1 == k2) return _listv[index(_m, i1, j1, k1)];

    // otherwise join the lists for all the cells in the range
    std::set<const SmoothedParticle*> joined;
    for (int i = i1; i <= i2; i++)
        for (int j = j1; j <= j2; j++)
            for (int k = k1; k <= k2; k++)
            {
                // add the particles for this list to the result set (if not already present)
                const vector<const SmoothedParticle*>& particles = _listv[index(_m, i, j, k)];
                joined.insert(particles.begin(), particles.end());
            }

    // copy the result back into a vector
    return vector<const SmoothedParticle*>(joined.begin(), joined.end());
}

////////////////////////////////////////////////////////////////////

const SmoothedParticle* SmoothedParticleGrid::nearestParticle(Vec r) const
{
    const SmoothedParticle* nearestParticle = nullptr;
    double nearestSquaredDistance = std::numeric_limits<double>::infinity();
    for (const SmoothedParticle* p : particlesFor(r))
    {
        double d2 = (r - p->center()).norm2();
        if (d2 < nearestSquaredDistance)
        {
            nearestParticle = p;
            nearestSquaredDistance = d2;
        }
    }
    return nearestParticle;
}

////////////////////////////////////////////////////////////////////
