/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "CellSnapshot.hpp"
#include "EntityCollection.hpp"
#include "Log.hpp"
#include "NR.hpp"
#include "Random.hpp"
#include "StringUtils.hpp"
#include "TextInFile.hpp"
#include "Units.hpp"

////////////////////////////////////////////////////////////////////

namespace
{
    // returns a Box object corresponding to the six elements in the given properties array starting at the given index
    Box box(const Array& prop, int i)
    {
        return Box(prop[i], prop[i + 1], prop[i + 2], prop[i + 3], prop[i + 4], prop[i + 5]);
    }

    // builds a smart grid in the specified spatial direction (box+0=x, box+1=y, box+2=z) and with the specified size,
    // and stores it in output parameter "grid", along with the minimum and maximum coordinate enclosing the cells
    void makegrid(const vector<Array>& propv, int dir, int gridsize, Array& grid, double& cmin, double& cmax)
    {
        int n = propv.size();

        // find the spatial range of the cells in the specified direction
        cmin = +std::numeric_limits<double>::infinity();
        cmax = -std::numeric_limits<double>::infinity();
        for (const Array& prop : propv)
        {
            cmin = min(cmin, prop[dir]);
            cmax = max(cmax, prop[dir + 3]);
        }

        // determine the cell distribution by binning at a decent resolution
        int nbins = gridsize * 100;
        double binwidth = (cmax - cmin) / nbins;
        vector<int> bins(nbins);
        for (const Array& prop : propv)
        {
            double center = 0.5 * (prop[dir] + prop[dir + 3]);
            bins[static_cast<int>((center - cmin) / binwidth)] += 1;
        }

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

    // returns the linear index for element (i,j,k) in a p*p*p table
    inline int index(int p, int i, int j, int k) { return ((i * p) + j) * p + k; }
}

////////////////////////////////////////////////////////////////////

// This is a helper class for organizing cuboidal cells in a smart grid, so that
// it is easy to retrieve the first cell that overlaps a given point in space.
// The Box object on which this class is based specifies a cuboid guaranteed to
// enclose all cells in the grid.
class CellSnapshot::CellGrid : public Box
{
    // data members initialized during construction
    const vector<Array>& _propv;   // reference to the original list of cells
    int _i;                        // the box index in the properties list of each cell
    int _p;                        // number of grid cells in each spatial direction
    Array _xgrid, _ygrid, _zgrid;  // the m+1 grid separation points for each spatial direction
    vector<vector<int>> _listv;    // the m*m*m lists of indices for cells overlapping each grid cell
    int _pmin, _pmax, _ptotal;     // minimum, maximum nr of cells in list; total nr of cells in listv

public:
    // The constructor creates a cuboidal grid of the specified number of grid cells in each
    // spatial direction, and for each of the grid cells it builds a list of all cells
    // overlapping the grid cell. In an attempt to distribute the cells evenly over the
    // grid cells, the sizes of the grid cells in each spatial direction are chosen so that
    // the cell centers are evenly distributed over the grid cells.
    CellGrid(const vector<Array>& propv, int boxindex, int gridsize) : _propv(propv), _i(boxindex), _p(gridsize)
    {
        // build the grids in each spatial direction
        double xmin, ymin, zmin, xmax, ymax, zmax;
        makegrid(propv, boxindex + 0, gridsize, _xgrid, xmin, xmax);
        makegrid(propv, boxindex + 1, gridsize, _ygrid, ymin, ymax);
        makegrid(propv, boxindex + 2, gridsize, _zgrid, zmin, zmax);
        setExtent(xmin, ymin, zmin, xmax, ymax, zmax);

        // make room for p*p*p grid cells
        _listv.resize(gridsize * gridsize * gridsize);

        // add each cell to the list for every grid cell that it overlaps
        int n = propv.size();
        for (int m = 0; m != n; ++m)
        {
            Box cell = box(propv[m], boxindex);

            // find indices for first and last grid cell overlapped by cell, in each spatial direction
            int i1 = NR::locateClip(_xgrid, cell.xmin());
            int i2 = NR::locateClip(_xgrid, cell.xmax());
            int j1 = NR::locateClip(_ygrid, cell.ymin());
            int j2 = NR::locateClip(_ygrid, cell.ymax());
            int k1 = NR::locateClip(_zgrid, cell.zmin());
            int k2 = NR::locateClip(_zgrid, cell.zmax());

            // add the cell to all grid cells in that 3D range
            for (int i = i1; i <= i2; i++)
                for (int j = j1; j <= j2; j++)
                    for (int k = k1; k <= k2; k++)
                    {
                        _listv[index(gridsize, i, j, k)].push_back(m);
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

    // This function returns the smallest number of cells overlapping a single grid cell.
    int minCellRefsPerCell() const { return _pmin; }

    // This function returns the largest number of cells overlapping a single grid cell.
    int maxCellRefsPerCell() const { return _pmax; }

    // This function returns the total number of cell references for all cells in the grid.
    int totalCellRefs() const { return _ptotal; }

    // This function returns the index (in the list originally passed to the constructor)
    // of the first cell in the list that overlaps the specified position,
    // or -1 if none of the cells in the list overlap the specified position.
    int cellIndexFor(Position r) const
    {
        // locate the grid cell containing the specified position
        int i = NR::locateClip(_xgrid, r.x());
        int j = NR::locateClip(_ygrid, r.y());
        int k = NR::locateClip(_zgrid, r.z());

        // search the list of cells for that grid cell
        for (int m : _listv[index(_p, i, j, k)])
        {
            if (box(_propv[m], _i).contains(r)) return m;
        }
        return -1;
    }

    // This function replaces the contents of the specified entity collection by the set of cells
    // that overlap the path with specified starting point and direction.
    // The weight of a cell is given by the length of the path segment inside the cell.
    void getEntities(EntityCollection& entities, Position bfr, Direction bfk) const
    {
        // use the values in these variables only after an intersection call that returns true
        double smin, smax;

        // initialize the output collection
        entities.clear();

        // verify that the path intersects the domain
        if (extent().intersects(bfr, bfk, smin, smax))
        {
            // find the indices for first and last grid cell, in each spatial direction,
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

            // loop over all grid cells in that 3D range
            for (int i = i1; i <= i2; i++)
                for (int j = j1; j <= j2; j++)
                    for (int k = k1; k <= k2; k++)
                    {
                        // if the path intersects the grid cell
                        Box gridcellbox(_xgrid[i], _ygrid[j], _zgrid[k], _xgrid[i + 1], _ygrid[j + 1], _zgrid[k + 1]);
                        if (gridcellbox.intersects(bfr, bfk, smin, smax))
                        {
                            // loop over all cells in this grid cell
                            for (int m : _listv[index(_p, i, j, k)])
                            {
                                // if the path intersects the cell, add the cell to the output collection
                                if (box(_propv[m], _i).intersects(bfr, bfk, smin, smax)) entities.add(m, smax - smin);
                            }
                        }
                    }
        }
    }
};

////////////////////////////////////////////////////////////////////

CellSnapshot::~CellSnapshot()
{
    delete _grid;
}

////////////////////////////////////////////////////////////////////

void CellSnapshot::readAndClose()
{
    // read the snapshot cell info into memory
    _propv = infile()->readAllRows();

    // close the file
    Snapshot::readAndClose();

    // inform the user
    log()->info("  Number of cells: " + std::to_string(_propv.size()));

    // if a mass density policy has been set, calculate masses and densities for all cells
    if (hasMassDensityPolicy())
    {
        // allocate vectors for mass and density
        size_t n = _propv.size();
        Array Mv(n);
        _rhov.resize(n);

        // get the maximum temperature, or zero of there is none
        double maxT = useTemperatureCutoff() ? maxTemperature() : 0.;

        // initialize statistics
        double totalOriginalMass = 0;
        double totalMetallicMass = 0;
        double totalEffectiveMass = 0;

        // loop over all sites/cells
        int numIgnored = 0;
        for (size_t m = 0; m != n; ++m)
        {
            const Array& prop = _propv[m];

            // original mass is zero if temperature is above cutoff or if imported mass/density is not positive
            double originalMass = 0.;
            if (maxT && prop[temperatureIndex()] > maxT)
                numIgnored++;
            else
                originalMass = max(0., massIndex() >= 0 ? prop[massIndex()]
                                                        : prop[densityIndex()] * box(prop, boxIndex()).volume());

            double metallicMass = originalMass * (useMetallicity() ? prop[metallicityIndex()] : 1.);
            double effectiveMass = metallicMass * multiplier();

            Mv[m] = effectiveMass;
            _rhov[m] = effectiveMass / box(prop, boxIndex()).volume();

            totalOriginalMass += originalMass;
            totalMetallicMass += metallicMass;
            totalEffectiveMass += effectiveMass;
        }

        // log mass statistics
        logMassStatistics(numIgnored, totalOriginalMass, totalMetallicMass, totalEffectiveMass);

        // remember the effective mass
        _mass = totalEffectiveMass;

        // construct a vector with the normalized cumulative cell densities
        if (n) NR::cdf(_cumrhov, Mv);
    }

    // if needed, construct a 3D-grid over the domain, and create a list of cells that overlap each grid cell
    if (hasMassDensityPolicy() || needGetEntities())
    {
        int gridsize = max(20, static_cast<int>(pow(_propv.size(), 1. / 3.) / 5));
        string size = std::to_string(gridsize);
        log()->info("Constructing intermediate " + size + "x" + size + "x" + size + " grid for cells...");
        _grid = new CellGrid(_propv, boxIndex(), gridsize);
        log()->info("  Smallest number of cells per grid cell: " + std::to_string(_grid->minCellRefsPerCell()));
        log()->info("  Largest  number of cells per grid cell: " + std::to_string(_grid->maxCellRefsPerCell()));
        log()->info("  Average  number of cells per grid cell: "
                    + StringUtils::toString(_grid->totalCellRefs() / double(gridsize * gridsize * gridsize), 'f', 1));
    }
}

////////////////////////////////////////////////////////////////////

Box CellSnapshot::extent() const
{
    // if there are no cells, return an empty box
    if (_propv.empty()) return Box();

    // if there is a cell grid, ask it to return the extent (it is already calculated)
    if (_grid) return _grid->extent();

    // otherwise find the spatial range of the cells
    double xmin = +std::numeric_limits<double>::infinity();
    double xmax = -std::numeric_limits<double>::infinity();
    double ymin = +std::numeric_limits<double>::infinity();
    double ymax = -std::numeric_limits<double>::infinity();
    double zmin = +std::numeric_limits<double>::infinity();
    double zmax = -std::numeric_limits<double>::infinity();
    for (const Array& prop : _propv)
    {
        xmin = min(xmin, prop[boxIndex() + 0]);
        xmax = max(xmax, prop[boxIndex() + 3]);
        ymin = min(ymin, prop[boxIndex() + 1]);
        ymax = max(ymax, prop[boxIndex() + 4]);
        zmin = min(zmin, prop[boxIndex() + 2]);
        zmax = max(zmax, prop[boxIndex() + 5]);
    }
    return Box(xmin, ymin, zmin, xmax, ymax, zmax);
}

////////////////////////////////////////////////////////////////////

int CellSnapshot::numEntities() const
{
    return _propv.size();
}

////////////////////////////////////////////////////////////////////

double CellSnapshot::volume(int m) const
{
    return box(_propv[m], boxIndex()).volume();
}

////////////////////////////////////////////////////////////////////

double CellSnapshot::density(int m) const
{
    return _rhov[m];
}

////////////////////////////////////////////////////////////////////

double CellSnapshot::density(Position bfr) const
{
    int m = _grid ? _grid->cellIndexFor(bfr) : -1;
    return m >= 0 ? _rhov[m] : 0.;
}

////////////////////////////////////////////////////////////////////

double CellSnapshot::mass() const
{
    return _mass;
}

////////////////////////////////////////////////////////////////////

Position CellSnapshot::position(int m) const
{
    return Position(box(_propv[m], boxIndex()).center());
}

////////////////////////////////////////////////////////////////////

Position CellSnapshot::generatePosition(int m) const
{
    return random()->position(box(_propv[m], boxIndex()));
}

////////////////////////////////////////////////////////////////////

Position CellSnapshot::generatePosition() const
{
    // if there are no cells, return the origin
    if (_propv.empty()) return Position();

    // select a cell according to its mass contribution
    int m = NR::locateClip(_cumrhov, random()->uniform());

    return generatePosition(m);
}

////////////////////////////////////////////////////////////////////

const Array& CellSnapshot::properties(int m) const
{
    return _propv[m];
}

////////////////////////////////////////////////////////////////////

int CellSnapshot::nearestEntity(Position bfr) const
{
    return _grid ? _grid->cellIndexFor(bfr) : -1;
}

////////////////////////////////////////////////////////////////////

void CellSnapshot::getEntities(EntityCollection& entities, Position bfr) const
{
    entities.addSingle(_grid->cellIndexFor(bfr));
}

////////////////////////////////////////////////////////////////////

void CellSnapshot::getEntities(EntityCollection& entities, Position bfr, Direction bfk) const
{
    _grid->getEntities(entities, bfr, bfk);
}

////////////////////////////////////////////////////////////////////
