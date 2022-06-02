/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "ParticleSnapshot.hpp"
#include "EntityCollection.hpp"
#include "Log.hpp"
#include "NR.hpp"
#include "Random.hpp"
#include "SmoothingKernel.hpp"
#include "StringUtils.hpp"
#include "TextInFile.hpp"
#include "Units.hpp"

////////////////////////////////////////////////////////////////////

/** Particle is an almost trivial helper class for working with smoothed particles, usually
    imported from a smoothed particle hydrodynamical (SPH) simulation. A Particle object
    holds the particle's center position, smoothing length, and effective mass (after any
    multipliers have been applied). Isolating these properties in a simple object allows efficient
    storage and retrieval, for example when organizing smoothed particles in a search grid. */
class ParticleSnapshot::Particle
{
public:
    /** The constructor arguments specify the particle attributes: particle index, coordinates of
        the center, smoothing length, and effective mass. */
    Particle(int m, double x, double y, double z, double h, double M) : _r{x, y, z}, _h(h), _M(M), _m(m) {}

    /** This function returns the index of the particle. */
    int index() const { return _m; }

    /** This function returns the coordinates of the center of the particle. */
    Vec center() const { return Vec(_r[0], _r[1], _r[2]); }

    /** This function returns the x, y, or z-coordinate of the center of the particle, depending on
        the value of \em dir: 1->x, 2->y, 3->z. For any other value of \em dir the behavior is
        undefined. */
    double center(int dir) const { return _r[dir - 1]; }

    /** This function returns the smoothing length of the particle. */
    double radius() const { return _h; }

    /** This function returns the mass of the particle. */
    double mass() const { return _M; }

    /** This function returns the effective volume of the particle. */
    double volume() const { return _h * _h * _h; }

    /** This function returns the effective density of the particle. */
    double density() const { return _M / (_h * _h * _h); }

    /** This function returns the impact radius of the path specified by a starting position and
        direction relative to the particle center (i.e. the minimum distance between the path and
        the center). If the path's starting position is located beyond the impact point (the point
        of minimum distance), the function returns infinity to indicate that the path does not
        intersect the particle. */
    double impact(Vec r, Vec k) const
    {
        // assume that k is normalized so we don't need to divide by its norm
        double s = Vec::dot(k, center() - r);
        if (s < 0.) return std::numeric_limits<double>::infinity();
        return (center() - (r + s * k)).norm();
    }

private:
    // data members received as constructor arguments
    double _r[3];  // center coordinates
    double _h;     // smoothing length
    double _M;     // total mass
    int _m;        // index
};

////////////////////////////////////////////////////////////////////

/** ParticleGrid is a helper class for organizing Particle instances in a smart
    grid, so that it is easy to retrieve a list of all particles that may overlap a particular
    point in space. The Box object on which this class is based specifies a cuboid guaranteed to
    enclose all particles in the grid. */
class ParticleSnapshot::ParticleGrid : public Box
{
private:
    // returns the linear index for cell (i,j,k) in a m*m*m table
    static inline int index(int m, int i, int j, int k) { return ((i * m) + j) * m + k; }

    // builds a smart grid in the specified spatial direction (1=x, 2=y, 3=z) and with the specified size,
    // and stores it in output parameter "grid", along with the minimum and maximum coordinate enclosing the particles
    static void makegrid(const vector<Particle>& pv, int dir, int gridsize, Array& grid, double& cmin, double& cmax)
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

public:
    /** The constructor creates a cuboidal grid of the specified number of grid cells in each
        spatial direction, and for each of the cells it builds a list of all particles (partially
        or fully) overlapping the cell. In an attempt to distribute the particles evenly over the
        cells, the sizes of the grid cells in each spatial direction are chosen so that the
        particle centers are evenly distributed over the cells. The internal particle lists store
        pointers to the particle objects contained in the provided list \em pv, so that list
        must not be modified or deallocated as long as this grid instance exists. */
    ParticleGrid(const vector<Particle>& pv, int gridsize) : _m(gridsize)
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
                        Box cell(_xgrid[i], _ygrid[j], _zgrid[k], _xgrid[i + 1], _ygrid[j + 1], _zgrid[k + 1]);
                        if (cell.intersects(rc, h)) _listv[index(gridsize, i, j, k)].push_back(&pv[p]);
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

    /** This function returns the smallest number of particles overlapping a single cell. */
    int minParticlesPerCell() const { return _pmin; }

    /** This function returns the largest number of particles overlapping a single cell. */
    int maxParticlesPerCell() const { return _pmax; }

    /** This function returns the total number of particle references for all cells in the grid. */
    int totalParticles() const { return _ptotal; }

    /** This function returns a list of all particles that may overlap the specified position. It
        locates the cell containing the specified position and returns the list of particles
        overlapping that cell. Thus the list may include particles that don't actually overlap the
        specified position. The objective of this class is to make this function very fast, while
        limiting the number of unnecessary particles in the returned list. */
    const vector<const Particle*>& particlesFor(Vec r) const
    {
        int i = NR::locateClip(_xgrid, r.x());
        int j = NR::locateClip(_ygrid, r.y());
        int k = NR::locateClip(_zgrid, r.z());
        return _listv[index(_m, i, j, k)];
    }

    /** This function returns a pointer to the particle centered nearest to the specified position,
        or the null pointer if the specified position is outside of the grid. */
    const Particle* nearestParticle(Vec r) const
    {
        const Particle* nearestParticle = nullptr;
        double nearestSquaredDistance = std::numeric_limits<double>::infinity();
        for (const Particle* particle : particlesFor(r))
        {
            double d2 = (r - particle->center()).norm2();
            if (d2 < nearestSquaredDistance)
            {
                nearestParticle = particle;
                nearestSquaredDistance = d2;
            }
        }
        return nearestParticle;
    }

    /** This function replaces the contents of the specified entity collection by the set of
        particles with a smoothing kernel that overlaps the specified point \f${\bf{r}}\f$. The
        weight corresponding to each particle is set to the particle's smoothing kernel value at
        the given point. If the given point does not overlap any particle, the collection will be
        empty. */
    void getEntities(EntityCollection& entities, Vec bfr, const SmoothingKernel* kernel) const
    {
        entities.clear();
        for (const Particle* particle : particlesFor(bfr))
        {
            double h = particle->radius();
            double u = (bfr - particle->center()).norm() / h;
            // if the point is inside the particle, add the particle to the output collection
            if (u <= 1.)
            {
                double w = kernel->density(u);
                entities.add(particle->index(), w);
            }
        }
    }

    /** This function replaces the contents of the specified entity collection by the set of
        particles with a smoothing kernel that overlaps the specified path with starting point
        \f${\bf{r}}\f$ and direction \f${\bf{k}}\f$. The weight of each particle is given by the
        effective length seen by the path as it crosses the particle's smoothing kernel. If the
        path does not overlap any particle, the collection will be empty. */
    void getEntities(EntityCollection& entities, Vec bfr, Vec bfk, const SmoothingKernel* kernel) const
    {
        // use the values in these variables only after a box/path intersection call that returns true
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
                            // loop over all particles in this grid cell
                            for (const Particle* particle : _listv[index(_m, i, j, k)])
                            {
                                double h = particle->radius();
                                double q = particle->impact(bfr, bfk) / h;

                                // if the path intersects the particle, add the particle to the output collection
                                if (q < 1.)
                                {
                                    double w = kernel->columnDensity(q) * h;
                                    entities.add(particle->index(), w);
                                }
                            }
                        }
                    }
        }
    }

private:
    int _m;                                  // number of grid cells in each spatial direction
    Array _xgrid, _ygrid, _zgrid;            // the m+1 grid separation points for each spatial direction
    vector<vector<const Particle*>> _listv;  // the m*m*m lists of particles overlapping each grid cell
    int _pmin, _pmax, _ptotal;               // minimum, maximum nr of particles in list; total nr of particles in listv
};

////////////////////////////////////////////////////////////////////

ParticleSnapshot::ParticleSnapshot() {}

////////////////////////////////////////////////////////////////////

ParticleSnapshot::~ParticleSnapshot()
{
    delete _grid;
}

////////////////////////////////////////////////////////////////////

void ParticleSnapshot::readAndClose()
{
    // read the particle info into memory
    // if the user configured a temperature cutoff, we skip high-temperature particles
    // if the user configured a mass-density policy, we skip zero-mass particles
    int numTempIgnored = 0;
    int numMassIgnored = 0;
    Array row;
    while (infile()->readRow(row))
    {
        if (useTemperatureCutoff() && row[temperatureIndex()] > maxTemperature())
            numTempIgnored++;
        else if (hasMassDensityPolicy() && row[massIndex()] == 0)
            numMassIgnored++;
        else
            _propv.push_back(row);
    }

    // close the file
    Snapshot::readAndClose();

    // log the number of particles
    if (!numTempIgnored && !numMassIgnored)
    {
        log()->info("  Number of particles: " + std::to_string(_propv.size()));
    }
    else
    {
        if (numTempIgnored)
            log()->info("  Number of high-temperature particles ignored: " + std::to_string(numTempIgnored));
        if (numMassIgnored) log()->info("  Number of zero-mass particles ignored: " + std::to_string(numMassIgnored));
        log()->info("  Number of particles retained: " + std::to_string(_propv.size()));
    }

    // if a mass density policy has been set, calculate masses and densities for all cells
    if (hasMassDensityPolicy())
    {
        // build a list of compact smoothed particle objects that we can organize in a grid
        double totalOriginalMass = 0;
        double totalMetallicMass = 0;
        double totalEffectiveMass = 0;
        int numParticles = _propv.size();
        _pv.reserve(numParticles);
        for (int m = 0; m != numParticles; ++m)
        {
            const Array& prop = _propv[m];

            double originalMass = prop[massIndex()];
            double metallicMass = originalMass * (useMetallicity() ? prop[metallicityIndex()] : 1.);
            double effectiveMass = metallicMass * multiplier();

            _pv.emplace_back(m, prop[positionIndex() + 0], prop[positionIndex() + 1], prop[positionIndex() + 2],
                             prop[sizeIndex()], effectiveMass);

            totalOriginalMass += originalMass;
            totalMetallicMass += metallicMass;
            totalEffectiveMass += effectiveMass;
        }

        // log mass statistics
        logMassStatistics(0, totalOriginalMass, totalMetallicMass, totalEffectiveMass);

        // if one of the total masses is negative, suppress the complete mass distribution
        if (totalOriginalMass < 0 || totalMetallicMass < 0 || totalEffectiveMass < 0)
        {
            log()->warning("  Total imported mass is negative; suppressing the complete mass distribution");
            _propv.clear();
            _pv.clear();
            totalEffectiveMass = 0;
        }

        // remember the effective mass
        _mass = totalEffectiveMass;

        // construct a vector with the normalized cumulative particle densities
        if (!_pv.empty()) NR::cdf(_cumrhov, _pv.size(), [this](int i) { return _pv[i].mass(); });
    }

    // if needed, build a list of compact particles without masses that we can organize in a grid
    if (!hasMassDensityPolicy() && needGetEntities())
    {
        int numParticles = _propv.size();
        _pv.reserve(numParticles);
        for (int m = 0; m != numParticles; ++m)
        {
            const Array& prop = _propv[m];
            _pv.emplace_back(m, prop[positionIndex() + 0], prop[positionIndex() + 1], prop[positionIndex() + 2],
                             prop[sizeIndex()], 0.);
        }
    }

    // if needed, construct a 3D-grid over the domain, and create a list of particles that overlap each grid cell
    if (hasMassDensityPolicy() || needGetEntities())
    {
        int gridsize = max(20, static_cast<int>(pow(_pv.size(), 1. / 3.) / 5));
        string size = std::to_string(gridsize);
        log()->info("Constructing intermediate " + size + "x" + size + "x" + size + " grid for particles...");
        _grid = new ParticleGrid(_pv, gridsize);
        log()->info("  Smallest number of particles per cell: " + std::to_string(_grid->minParticlesPerCell()));
        log()->info("  Largest  number of particles per cell: " + std::to_string(_grid->maxParticlesPerCell()));
        log()->info("  Average  number of particles per cell: "
                    + StringUtils::toString(_grid->totalParticles() / double(gridsize * gridsize * gridsize), 'f', 1));
    }
}

////////////////////////////////////////////////////////////////////

void ParticleSnapshot::setSmoothingKernel(const SmoothingKernel* kernel)
{
    _kernel = kernel;
}

////////////////////////////////////////////////////////////////////

Box ParticleSnapshot::extent() const
{
    // if there are no particles, return an empty box
    if (_propv.empty()) return Box();

    // if there is a particle grid, ask it to return the extent (it is already calculated)
    if (_grid) return _grid->extent();

    // otherwise find the spatial range of the particles assuming a finite support kernel
    double xmin = +std::numeric_limits<double>::infinity();
    double xmax = -std::numeric_limits<double>::infinity();
    double ymin = +std::numeric_limits<double>::infinity();
    double ymax = -std::numeric_limits<double>::infinity();
    double zmin = +std::numeric_limits<double>::infinity();
    double zmax = -std::numeric_limits<double>::infinity();
    for (const Array& prop : _propv)
    {
        xmin = min(xmin, prop[positionIndex() + 0] - prop[sizeIndex()]);
        xmax = max(xmax, prop[positionIndex() + 0] + prop[sizeIndex()]);
        ymin = min(ymin, prop[positionIndex() + 1] - prop[sizeIndex()]);
        ymax = max(ymax, prop[positionIndex() + 1] + prop[sizeIndex()]);
        zmin = min(zmin, prop[positionIndex() + 2] - prop[sizeIndex()]);
        zmax = max(zmax, prop[positionIndex() + 2] + prop[sizeIndex()]);
    }
    return Box(xmin, ymin, zmin, xmax, ymax, zmax);
}

////////////////////////////////////////////////////////////////////

int ParticleSnapshot::numEntities() const
{
    return _propv.size();
}

////////////////////////////////////////////////////////////////////

double ParticleSnapshot::volume(int m) const
{
    return _pv[m].volume();
}

////////////////////////////////////////////////////////////////////

double ParticleSnapshot::density(int m) const
{
    return _pv[m].density();
}

////////////////////////////////////////////////////////////////////

double ParticleSnapshot::density(Position bfr) const
{
    double sum = 0.;
    if (_grid)
        for (const Particle* p : _grid->particlesFor(bfr))
        {
            double h = p->radius();
            double u = (bfr - p->center()).norm() / h;
            sum += _kernel->density(u) * p->density();
        }
    return sum > 0. ? sum : 0.;  // guard against negative densities
}

////////////////////////////////////////////////////////////////////

double ParticleSnapshot::mass() const
{
    return _mass;
}

////////////////////////////////////////////////////////////////////

Position ParticleSnapshot::position(int m) const
{
    return Position(_propv[m][positionIndex() + 0], _propv[m][positionIndex() + 1], _propv[m][positionIndex() + 2]);
}

////////////////////////////////////////////////////////////////////

Position ParticleSnapshot::generatePosition(int m) const
{
    // get center position and size for this particle
    Position rc(_propv[m][positionIndex() + 0], _propv[m][positionIndex() + 1], _propv[m][positionIndex() + 2]);
    double h = _propv[m][sizeIndex()];

    // sample random position inside the smoothed unit volume
    double u = _kernel->generateRadius();
    Direction k = random()->direction();

    return Position(rc + k * u * h);
}

////////////////////////////////////////////////////////////////////

Position ParticleSnapshot::generatePosition() const
{
    // if there are no particles, return the origin
    if (_propv.empty()) return Position();

    // select a particle according to its mass contribution
    int m = NR::locateClip(_cumrhov, random()->uniform());

    return generatePosition(m);
}

////////////////////////////////////////////////////////////////////

const Array& ParticleSnapshot::properties(int m) const
{
    return _propv[m];
}

////////////////////////////////////////////////////////////////////

int ParticleSnapshot::nearestEntity(Position bfr) const
{
    const Particle* nearestParticle = _grid ? _grid->nearestParticle(bfr) : nullptr;
    return nearestParticle ? nearestParticle->index() : -1;
}

////////////////////////////////////////////////////////////////////

void ParticleSnapshot::getEntities(EntityCollection& entities, Position bfr) const
{
    _grid->getEntities(entities, bfr, _kernel);
}

////////////////////////////////////////////////////////////////////

void ParticleSnapshot::getEntities(EntityCollection& entities, Position bfr, Direction bfk) const
{
    _grid->getEntities(entities, bfr, bfk, _kernel);
}

////////////////////////////////////////////////////////////////////
