/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "TetraMeshSnapshot.hpp"
#include "EntityCollection.hpp"
#include "FatalError.hpp"
#include "Log.hpp"
#include "NR.hpp"
#include "Parallel.hpp"
#include "ParallelFactory.hpp"
#include "PathSegmentGenerator.hpp"
#include "ProcessManager.hpp"
#include "Random.hpp"
#include "SiteListInterface.hpp"
#include "SpatialGridPath.hpp"
#include "SpatialGridPlotFile.hpp"
#include "StringUtils.hpp"
#include "Table.hpp"
#include "TextInFile.hpp"
#include "Units.hpp"
#include <iostream>
#include <set>
#include "container.hh"

////////////////////////////////////////////////////////////////////

// class to hold the information about a Tetra cell that is relevant for calculating paths and densities
class TetraMeshSnapshot::Cell
{
public:
    Vec _r;                  // site position
    vector<int> _neighbors;  // list of neighbor indices in _cells vector
    Array _properties;       // user-defined properties, if any

public:
    // constructor stores the specified site position; the other data members are set to zero or empty
    Cell(Vec r) : _r(r) {}

    // constructor derives the site position from the first three property values and stores the user properties;
    // the other data members are set to zero or empty
    Cell(const Array& prop) : _r(prop[0], prop[1], prop[2]), _properties{prop} {}  // WIP

    // adjusts the site position with the specified offset
    void relax(double cx, double cy, double cz) { _r += Vec(cx, cy, cz); }

    // initializes the receiver with information taken from the specified fully computed Tetra cell
    void init(voro::voronoicell_neighbor& cell) { cell.neighbors(_neighbors); }

    // clears the information taken from a Tetra cell so it can be reinitialized
    void clear() { _neighbors.clear(); }

    // returns the cell's site position
    Vec position() const { return _r; }

    // returns the x coordinate of the cell's site position
    double x() const { return _r.x(); }

    // returns the squared distance from the cell's site to the specified point
    double squaredDistanceTo(Vec r) const { return (r - _r).norm2(); }

    // returns a list of neighboring cell/site ids
    const vector<int>& neighbors() { return _neighbors; }

    // returns the cell/site user properties, if any
    const Array& properties() { return _properties; }
};

////////////////////////////////////////////////////////////////////

TetraMeshSnapshot::Plucker::Plucker() {}

////////////////////////////////////////////////////////////////////

TetraMeshSnapshot::Plucker::Plucker(const Vec& pos, const Vec& dir) : U(dir), V(Vec::cross(dir, pos)) {}

////////////////////////////////////////////////////////////////////

// permuted inner product
inline double TetraMeshSnapshot::Plucker::dot(const Plucker& a, const Plucker& b)
{
    return Vec::dot(a.U, b.V) + Vec::dot(b.U, a.V);
}

////////////////////////////////////////////////////////////////////

TetraMeshSnapshot::Edge::Edge(const int i1, const int i2, const Vec& v1, const Vec& v2)
    : Plucker((i1 < i2) ? v1 : v2, (i1 < i2) ? v2 : v1), i1(std::min(i1, i2)), i2(std::max(i1, i2))
{}

////////////////////////////////////////////////////////////////////

bool TetraMeshSnapshot::Edge::operator==(const Edge& edge) const
{
    return (i1 == edge.i1 && i2 == edge.i2);
}

////////////////////////////////////////////////////////////////////

size_t TetraMeshSnapshot::Edge::HashFunction::operator()(const Edge* edge) const
{
    return ((size_t)edge->i2 << 32) + (size_t)edge->i1;
}

////////////////////////////////////////////////////////////////////

bool TetraMeshSnapshot::Edge::EqualFunction::operator()(const Edge* e1, const Edge* e2) const
{
    return *e1 == *e2;
}

////////////////////////////////////////////////////////////////////

TetraMeshSnapshot::Tetra::Tetra(const vector<Cell*>& _cells, int i, int j, int k, int l)
{
    _vertices[0] = &_cells[i]->_r;
    _vertices[1] = &_cells[j]->_r;
    _vertices[2] = &_cells[k]->_r;
    _vertices[3] = &_cells[l]->_r;
    _vertex_indices[0] = i;
    _vertex_indices[1] = j;
    _vertex_indices[2] = k;
    _vertex_indices[3] = l;

    double xmin = DBL_MAX;
    double ymin = DBL_MAX;
    double zmin = DBL_MAX;
    double xmax = -DBL_MAX;
    double ymax = -DBL_MAX;
    double zmax = -DBL_MAX;
    for (const Vec* vertex : _vertices)
    {
        xmin = min(xmin, vertex->x());
        ymin = min(ymin, vertex->y());
        zmin = min(zmin, vertex->z());
        xmax = max(xmax, vertex->x());
        ymax = max(ymax, vertex->y());
        zmax = max(zmax, vertex->z());
    }
    setExtent(Box(xmin, ymin, zmin, xmax, ymax, zmax));

    _volume = Box::volume();  // WIP
}

////////////////////////////////////////////////////////////////////

int TetraMeshSnapshot::Tetra::getEdge(int t1, int t2) const
{
    // negative index shows that the edge is inverted
    if (t1 > t2)
        return -((t2 == 0) ? t1 - 1 : t2 + t1);
    else
        return (t1 == 0) ? t2 - 1 : t1 + t2;
}

////////////////////////////////////////////////////////////////////

double TetraMeshSnapshot::Tetra::orient()
{
    const Vec e01 = *_vertices[1] - *_vertices[0];
    const Vec e02 = *_vertices[2] - *_vertices[0];
    const Vec e03 = *_vertices[3] - *_vertices[0];
    // this convention makes edges go clockwise around leaving rays from inside the tetrahedron
    // so their plucker products are all negative if the ray leaves
    double orientation = Vec::dot(Vec::cross(e01, e02), e03);
    if (orientation > 0)
    {
        // swap last 2, this means first 2 indices can be ordered i < j
        std::swap(_vertices[2], _vertices[3]);
        std::swap(_vertex_indices[2], _vertex_indices[3]);
    }
    return orientation;
}

////////////////////////////////////////////////////////////////////

void TetraMeshSnapshot::Tetra::addEdges(EdgeSet& edgeset)
{
    for (int t1 = 0; t1 < 3; t1++)
    {
        for (int t2 = 1; t2 < 4; t2++)
        {
            if (t2 <= t1) continue;

            Edge* edge = new Edge(_vertex_indices[t1], _vertex_indices[t2], *_vertices[t1], *_vertices[t2]);
            auto res = edgeset.insert(edge);

            // insertion unsuccessful (edge already exists)
            if (!res.second)
            {
                delete edge;
                edge = *res.first;
            }

            _edges[getEdge(t1, t2)] = edge;
        }
    }
}

////////////////////////////////////////////////////////////////////

int TetraMeshSnapshot::Tetra::shareFace(const Tetra* other) const
{
    int equal_vert = 0;
    int opposite = 0 + 1 + 2 + 3;

    for (int i = 0; i < 4; i++)
    {
        int vi = _vertex_indices[i];

        for (int j = 0; j < 4; j++)
        {
            int vj = other->_vertex_indices[j];

            if (vi == vj)
            {
                equal_vert++;
                opposite -= i;
                break;
            }
        }
    }
    if (equal_vert == 3) return opposite;
    return -1;
}

////////////////////////////////////////////////////////////////////

bool TetraMeshSnapshot::Tetra::intersects(std::array<double, 3>& prods, const Plucker& ray, int face,
                                          bool leaving = true) const
{
    int t0 = (face + 1) % 4;
    int t1 = (face + 2) % 4;
    int t2 = (face + 3) % 4;

    // if face is even we should swap two edges
    if (face % 2 == 0) std::swap(t0, t2);

    Edge* edges[3] = {_edges[getEdge(t1, t2)], _edges[getEdge(t2, t0)], _edges[getEdge(t0, t1)]};

    double sum = 0;
    for (int i = 0; i < 3; i++)
    {
        double prod = Plucker::dot(ray, *_edges[i]);
        if (leaving != (prod < 0)) return false;
        prods[i] = prod;
        sum += prod;
    }
    prods = {prods[0] / sum, prods[1] / sum, prods[2] / sum};
    return true;
}

////////////////////////////////////////////////////////////////////

bool TetraMeshSnapshot::Tetra::equals(const Tetra* other) const
{
    for (int v : _vertex_indices)
    {
        bool match = false;
        for (int u : other->_vertex_indices)
        {
            if (u == v)
            {
                match = true;
                break;
            }
        }
        if (!match) return false;
    }
    return true;
}

////////////////////////////////////////////////////////////////////

bool TetraMeshSnapshot::Tetra::SameSide(const Vec& v0, const Vec& v1, const Vec& v2, const Vec& v3,
                                        const Vec& pos) const
{
    Vec normal = Vec::cross(v1 - v0, v2 - v0);
    double dotV4 = Vec::dot(normal, v3 - v0);
    double dotP = Vec::dot(normal, pos - v0);
    return (dotV4 > 0) == (dotP > 0);
}

////////////////////////////////////////////////////////////////////

bool TetraMeshSnapshot::Tetra::inside(const Position& bfr) const
{
    // very poor implementation
    const Vec& v0 = *_vertices[0];
    const Vec& v1 = *_vertices[1];
    const Vec& v2 = *_vertices[2];
    const Vec& v3 = *_vertices[3];
    return SameSide(v0, v1, v2, v3, bfr) && SameSide(v1, v2, v3, v0, bfr) && SameSide(v2, v3, v0, v1, bfr)
           && SameSide(v3, v0, v1, v2, bfr);
}

////////////////////////////////////////////////////////////////////

Position TetraMeshSnapshot::Tetra::position() const
{
    Position pos;
    for (int i = 0; i < 4; i++) pos += *_vertices[i];
    pos /= 4;
    return pos;
}

////////////////////////////////////////////////////////////////////

double TetraMeshSnapshot::Tetra::volume() const
{
    return _volume;
}

////////////////////////////////////////////////////////////////////

const Array& TetraMeshSnapshot::Tetra::properties()
{
    return _properties;
}

////////////////////////////////////////////////////////////////////

namespace
{
    template<typename T> bool invec(const vector<T>& vec, const T& e)
    {
        return std::find(vec.begin(), vec.end(), e) != vec.end();
    }

    double det4(double x0, double y0, double z0, double x1, double y1, double z1, double x2, double y2, double z2,
                double x3, double y3, double z3)
    {
        return x0 * y1 * z2 - x0 * y1 * z3 - x0 * y2 * z1 + x0 * y2 * z3 + x0 * y3 * z1 - x0 * y3 * z2 - x1 * y0 * z2
               + x1 * y0 * z3 + x1 * y2 * z0 - x1 * y2 * z3 - x1 * y3 * z0 + x1 * y3 * z2 + x2 * y0 * z1 - x2 * y0 * z3
               - x2 * y1 * z0 + x2 * y1 * z3 + x2 * y3 * z0 - x2 * y3 * z1 - x3 * y0 * z1 + x3 * y0 * z2 + x3 * y1 * z0
               - x3 * y1 * z2 - x3 * y2 * z0 + x3 * y2 * z1;
    }
}

////////////////////////////////////////////////////////////////////

TetraMeshSnapshot::TetraMeshSnapshot() {}

////////////////////////////////////////////////////////////////////

TetraMeshSnapshot::~TetraMeshSnapshot()
{
    for (auto cell : _cells) delete cell;
    for (auto tetra : _tetrahedra) delete tetra;
    for (auto tree : _blocktrees) delete tree;
}

////////////////////////////////////////////////////////////////////

void TetraMeshSnapshot::setExtent(const Box& extent)
{
    _extent = extent;
    _eps = 1e-12 * extent.widths().norm();
}

////////////////////////////////////////////////////////////////////

TetraMeshSnapshot::TetraMeshSnapshot(const SimulationItem* item, const Box& extent, string filename, bool relax)
{
    // read the input file
    TextInFile in(item, filename, "Tetra sites");
    in.addColumn("position x", "length", "pc");
    in.addColumn("position y", "length", "pc");
    in.addColumn("position z", "length", "pc");
    Array coords;
    while (in.readRow(coords)) _cells.push_back(new Cell(Vec(coords[0], coords[1], coords[2])));
    in.close();

    // calculate the Tetra cells
    setContext(item);
    setExtent(extent);
    buildMesh(relax);
    // buildSearchPerBlock();
}

////////////////////////////////////////////////////////////////////

TetraMeshSnapshot::TetraMeshSnapshot(const SimulationItem* item, const Box& extent, SiteListInterface* sli, bool relax)
{
    // prepare the data
    int n = sli->numSites();
    _cells.resize(n);
    for (int m = 0; m != n; ++m) _cells[m] = new Cell(sli->sitePosition(m));

    // calculate the Tetra cells
    setContext(item);
    setExtent(extent);
    buildMesh(relax);
    // buildSearchPerBlock();
}

////////////////////////////////////////////////////////////////////

TetraMeshSnapshot::TetraMeshSnapshot(const SimulationItem* item, const Box& extent, const vector<Vec>& sites,
                                     bool relax)
{
    // prepare the data
    int n = sites.size();
    _cells.resize(n);
    for (int m = 0; m != n; ++m) _cells[m] = new Cell(sites[m]);

    // calculate the Tetra cells
    setContext(item);
    setExtent(extent);
    buildMesh(relax);
    // buildSearchPerBlock();
}

////////////////////////////////////////////////////////////////////

namespace
{
    // maximum number of Tetra sites processed between two invocations of infoIfElapsed()
    const int logProgressChunkSize = 1000;

    // maximum number of Tetra grid construction iterations
    const int maxConstructionIterations = 5;

    // function to erase null pointers from a vector of pointers in one go; returns the new size
    template<class T> size_t eraseNullPointers(vector<T*>& v)
    {
        // scan to the first null (or to the end)
        auto from = v.begin();
        while (from != v.end() && *from) ++from;

        // copy the tail, overwriting any nulls
        auto to = from;
        while (from != v.end())
        {
            if (*from) *to++ = *from;
            ++from;
        }
        v.erase(to, v.end());
        return v.size();
    }
}

////////////////////////////////////////////////////////////////////

void TetraMeshSnapshot::buildMesh(bool relax)
{
    const double A = extent().rmax().norm() * 1000;  // WIP

    int numCells = _cells.size();

    // remove sites that lie outside of the domain
    int numOutside = 0;
    for (int m = 0; m != numCells; ++m)
    {
        if (!_extent.contains(_cells[m]->position()))
        {
            delete _cells[m];
            _cells[m] = 0;
            numOutside++;
        }
    }
    if (numOutside) numCells = eraseNullPointers(_cells);

    // sort sites in order of increasing x coordinate to accelerate search for nearby sites
    std::sort(_cells.begin(), _cells.end(), [](Cell* c1, Cell* c2) { return c1->x() < c2->x(); });

    // remove sites that lie too nearby another site
    int numNearby = 0;
    for (int m = 0; m != numCells; ++m)
    {
        for (int j = m + 1; j != numCells && _cells[j]->x() - _cells[m]->x() < _eps; ++j)
        {
            if ((_cells[j]->position() - _cells[m]->position()).norm2() < _eps * _eps)
            {
                delete _cells[m];
                _cells[m] = 0;
                numNearby++;
                break;
            }
        }
    }
    if (numNearby) numCells = eraseNullPointers(_cells);

    // log the number of sites
    if (!numOutside && !numNearby)
    {
        log()->info("  Number of sites: " + std::to_string(numCells));
    }
    else
    {
        if (numOutside) log()->info("  Number of sites outside domain: " + std::to_string(numOutside));
        if (numNearby) log()->info("  Number of sites too nearby others: " + std::to_string(numNearby));
        log()->info("  Number of sites retained: " + std::to_string(numCells));
    }

    // abort if there are no cells to calculate
    if (numCells <= 0) return;

    // calculate number of blocks in each direction based on number of cells
    _nb = max(3, min(250, static_cast<int>(cbrt(numCells))));
    _nb2 = _nb * _nb;
    _nb3 = _nb * _nb * _nb;

    // ========= RELAXATION =========

    // if requested, perform a single relaxation step
    if (relax)
    {
        // table to hold the calculate relaxation offset for each site
        // (initialized to zero so we can communicate the result between parallel processes using sumAll)
        Table<2> offsets(numCells, 3);

        // add the retained original sites to a temporary Voronoi container, using the cell index m as ID
        voro::container vcon(_extent.xmin(), _extent.xmax(), _extent.ymin(), _extent.ymax(), _extent.zmin(),
                             _extent.zmax(), _nb, _nb, _nb, false, false, false, 16);
        for (int m = 0; m != numCells; ++m)
        {
            Vec r = _cells[m]->position();
            vcon.put(m, r.x(), r.y(), r.z());
        }

        // compute the cell in the Voronoi tesselation corresponding to each site
        // and store the cell's centroid (relative to the site position) as the relaxation offset
        log()->info("Relaxing Voronoi tessellation with " + std::to_string(numCells) + " cells");
        log()->infoSetElapsed(numCells);
        auto parallel = log()->find<ParallelFactory>()->parallelDistributed();
        parallel->call(numCells, [this, &vcon, &offsets](size_t firstIndex, size_t numIndices) {
            // allocate a separate cell calculator for each thread to avoid conflicts
            voro::voro_compute<voro::container> vcompute(vcon, _nb, _nb, _nb);
            // allocate space for the resulting cell info
            voro::voronoicell vcell;

            // loop over all cells and work on the ones that have a particle index in our dedicated range
            // (we cannot access cells in the container based on cell index m without building an extra data structure)
            int numDone = 0;
            voro::c_loop_all vloop(vcon);
            if (vloop.start()) do
                {
                    size_t m = vloop.pid();
                    if (m >= firstIndex && m < firstIndex + numIndices)
                    {
                        // compute the cell and store its centroid as relaxation offset
                        bool ok = vcompute.compute_cell(vcell, vloop.ijk, vloop.q, vloop.i, vloop.j, vloop.k);
                        if (ok) vcell.centroid(offsets(m, 0), offsets(m, 1), offsets(m, 2));

                        // log message if the minimum time has elapsed
                        numDone = (numDone + 1) % logProgressChunkSize;
                        if (numDone == 0) log()->infoIfElapsed("Computed Voronoi cells: ", logProgressChunkSize);
                    }
                } while (vloop.inc());
            if (numDone > 0) log()->infoIfElapsed("Computed Voronoi cells: ", numDone);
        });

        // communicate the calculated offsets between parallel processes, if needed, and apply them to the cells
        ProcessManager::sumToAll(offsets.data());
        for (int m = 0; m != numCells; ++m) _cells[m]->relax(offsets(m, 0), offsets(m, 1), offsets(m, 2));
    }

    // ========= FINAL GRID =========

    // add the final sites to a temporary Tetra container, using the cell index m as ID
    voro::container vcon(-A, A, -A, A, -A, A, _nb, _nb, _nb, false, false, false, 16);
    for (int m = 0; m != numCells; ++m)
    {
        Vec r = _cells[m]->position();
        vcon.put(m, r.x(), r.y(), r.z());
    }

    // for each site:
    //   - compute the corresponding cell in the Tetra tesselation
    //   - extract and copy the relevant information to the cell object with the corresponding index in our vector
    log()->info("Constructing intermediate Voronoi tessellation with " + std::to_string(numCells) + " cells");
    log()->infoSetElapsed(numCells);
    // allocate a separate cell calculator for each thread to avoid conflicts
    voro::voro_compute<voro::container> vcompute(vcon, _nb, _nb, _nb);
    // allocate space for the resulting cell info
    voro::voronoicell_neighbor vcell;
    std::vector<voro::voronoicell_neighbor> vcells;

    // loop over all cells and work on the ones that have a particle index in our dedicated range
    // (we cannot access cells in the container based on cell index m without building an extra data structure)
    int numDone = 0;
    voro::c_loop_all vloop(vcon);
    if (vloop.start()) do
        {
            size_t m = vloop.pid();
            // compute the cell and copy all relevant information to the cell object that will stay around
            bool ok = vcompute.compute_cell(vcell, vloop.ijk, vloop.q, vloop.i, vloop.j, vloop.k);
            if (ok) _cells[m]->init(vcell);

            // log message if the minimum time has elapsed
            numDone = (numDone + 1) % logProgressChunkSize;
            if (numDone == 0) log()->infoIfElapsed("Computed Tetra cells: ", logProgressChunkSize);
        } while (vloop.inc());
    if (numDone > 0) log()->infoIfElapsed("Computed Tetra cells: ", numDone);

    for (int i = 0; i < numCells; i++)
    {
        Cell* c1 = _cells[i];

        for (int j : c1->neighbors())
        {
            // first 2 indices can always be ordered i < j
            // last 2 indices can be swapped if orientation is not correct
            if (j < 0 || j > i) continue;

            Cell* c2 = _cells[j];

            for (int k : c2->neighbors())
            {
                if (k < 0 || k == i) continue;

                // mutual neighbours
                if (!invec(c1->neighbors(), k)) continue;

                for (int l : _cells[k]->neighbors())
                {
                    if (l < 0 || l == j || l == i) continue;

                    // mutual neighbours of both c1 and c2
                    if (!invec(c1->neighbors(), l) || !invec(c2->neighbors(), l)) continue;

                    // no duplicates in different order // probably use set
                    Tetra* tetra = new Tetra(_cells, i, j, k, l);
                    if (inTetrahedra(tetra)) continue;

                    // check if tetrahedron is Delaunay
                    Vec v0 = _cells[i]->_r;
                    Vec v1 = _cells[j]->_r;
                    Vec v2 = _cells[k]->_r;
                    Vec v3 = _cells[l]->_r;
                    double x0 = v0.x(), y0 = v0.y(), z0 = v0.z();
                    double x1 = v1.x(), y1 = v1.y(), z1 = v1.z();
                    double x2 = v2.x(), y2 = v2.y(), z2 = v2.z();
                    double x3 = v3.x(), y3 = v3.y(), z3 = v3.z();

                    double r0 = v0.norm2();
                    double r1 = v1.norm2();
                    double r2 = v2.norm2();
                    double r3 = v3.norm2();

                    double a = det4(x0, y0, z0, x1, y1, z1, x2, y2, z2, x3, y3, z3);
                    double Dx = det4(r0, y0, z0, r1, y1, z1, r2, y2, z2, r3, y3, z3);
                    double Dy = -det4(r0, x0, z0, r1, x1, z1, r2, x2, z2, r3, x3, z3);
                    double Dz = det4(r0, x0, y0, r1, x1, y1, r2, x2, y2, r3, x3, y3);

                    Vec center(Dx / (2 * a), Dy / (2 * a), Dz / (2 * a));
                    double R = (center - v0).norm2();

                    bool delaunay = true;
                    for (int v : tetra->_vertex_indices)
                    {
                        Cell* c = _cells[v];
                        for (int n : c->neighbors())
                        {
                            if (n < 0 || n == i || n == j || n == k || n == l) continue;
                            double r = (center - _cells[n]->_r).norm2();
                            if (r < R)
                            {
                                delaunay = false;
                                goto end_delaunay;
                            }
                        }
                    }
                end_delaunay:
                    if (!delaunay) continue;

                    // orient tetrahedron in the same consistent way
                    tetra->orient();

                    //add edges
                    tetra->addEdges(_edgeset);

                    centers.push_back(center);
                    radii.push_back(sqrt(R));
                    _tetrahedra.push_back(tetra);
                }
            }
        }
    }

    // find neighbors brute force
    for (size_t i = 0; i < _tetrahedra.size(); i++)
    {
        Tetra* tetra = _tetrahedra[i];

        for (size_t j = 0; j < _tetrahedra.size(); j++)
        {
            if (i == j) continue;

            const Tetra* other = _tetrahedra[j];

            int shared = tetra->shareFace(other);
            if (shared != -1) tetra->_neighbors[shared] = j;
        }
    }
}

////////////////////////////////////////////////////////////////////

void TetraMeshSnapshot::calculateVolume()
{
    //WIP
}

////////////////////////////////////////////////////////////////////

void TetraMeshSnapshot::calculateDensityAndMass()
{
    // allocate vectors for mass and density
    int numCells = _tetrahedra.size();
    _rhov.resize(numCells);
    Array Mv(numCells);

    // get the maximum temperature, or zero of there is none
    double maxT = useTemperatureCutoff() ? maxTemperature() : 0.;

    // initialize statistics
    double totalOriginalMass = 0;
    double totalMetallicMass = 0;
    double totalEffectiveMass = 0;

    // loop over all sites/cells
    int numIgnored = 0;
    for (int m = 0; m != numCells; ++m)
    {
        const Array& prop = _tetrahedra[m]->properties();

        // original mass is zero if temperature is above cutoff or if imported mass/density is not positive
        double originalDensity = 0.;
        double originalMass = 0.;
        if (maxT && prop[temperatureIndex()] > maxT)
        {
            numIgnored++;
        }
        else
        {
            double volume = _tetrahedra[m]->volume();
            originalDensity = max(0., densityIndex() >= 0 ? prop[densityIndex()] : prop[massIndex()] / volume);
            originalMass = max(0., massIndex() >= 0 ? prop[massIndex()] : prop[densityIndex()] * volume);
        }

        double effectiveDensity = originalDensity * (useMetallicity() ? prop[metallicityIndex()] : 1.) * multiplier();
        double metallicMass = originalMass * (useMetallicity() ? prop[metallicityIndex()] : 1.);
        double effectiveMass = metallicMass * multiplier();

        _rhov[m] = effectiveDensity;
        Mv[m] = effectiveMass;

        totalOriginalMass += originalMass;
        totalMetallicMass += metallicMass;
        totalEffectiveMass += effectiveMass;
    }

    // log mass statistics
    logMassStatistics(numIgnored, totalOriginalMass, totalMetallicMass, totalEffectiveMass);

    // remember the effective mass
    _mass = totalEffectiveMass;

    // construct a vector with the normalized cumulative site densities
    if (numCells) NR::cdf(_cumrhov, Mv);
}

////////////////////////////////////////////////////////////////////

bool TetraMeshSnapshot::inTetrahedra(const Tetra* tetra) const
{
    for (const Tetra* t : _tetrahedra)
    {
        if (tetra->equals(t)) return true;
    }
    return false;
}

////////////////////////////////////////////////////////////////////

void TetraMeshSnapshot::writeGridPlotFiles(const SimulationItem* probe) const
{
    // create the plot files
    // SpatialGridPlotFile plotxy(probe, probe->itemName() + "_grid_xy");
    // SpatialGridPlotFile plotxz(probe, probe->itemName() + "_grid_xz");
    // SpatialGridPlotFile plotyz(probe, probe->itemName() + "_grid_yz");
    SpatialGridPlotFile plotxyz(probe, probe->itemName() + "_grid_xyz");

    // vector<double> test = {1, 0, 0, 0, 1, 0, 0, -1, 0, 0, 0, 1};
    // vector<int> test2 = {3, 0, 1, 2, 3, 0, 1, 3, 3, 0, 2, 3, 3, 1, 2, 3};
    // plotxyz.writePolyhedron(test, test2);
    // return;

    ////////////////////////////////////////////////////////////////////////////////////////////////
    std::ofstream outputFile("input.txt");
    outputFile << "voronois=" << _cells.size() << "\n";
    for (size_t i = 0; i < _cells.size(); i++)
    {
        outputFile << "voronoi=" << i << "\n";
        Vec& r = _cells[i]->_r;
        outputFile << r.x() << ", " << r.y() << ", " << r.z() << "\n";

        outputFile << _cells[i]->neighbors().size() << " neighbors=";
        for (int n : _cells[i]->neighbors())
        {
            if (n >= 0) outputFile << " " << n << ",";
        }
        outputFile << "\n";
    }

    outputFile << "tetrahedra=" << _tetrahedra.size() << "\n";
    for (size_t i = 0; i < _tetrahedra.size(); i++)
    {
        const Tetra* tetra = _tetrahedra[i];
        outputFile << "tetrahedron=" << i << "\nvertices=";
        for (size_t l = 0; l < 4; l++)
        {
            outputFile << " " << tetra->_vertex_indices[l] << ",";
        }

        outputFile << "\n" << tetra->_neighbors.size() << " neighbors=";
        for (size_t j = 0; j < 4; j++)
        {
            outputFile << " " << tetra->_neighbors[j] << ",";
        }
        outputFile << "\n";
    }
    for (size_t i = 0; i < _tetrahedra.size(); i++)
    {
        outputFile << "circumsphere=" << i << "\n";
        outputFile << centers[i].x() << "," << centers[i].y() << "," << centers[i].z() << "\n";
        outputFile << radii[i] << "\n";
    }

    // Close the output file
    outputFile.close();
    ////////////////////////////////////////////////////////////////////////////////////////////////

    for (size_t i = 0; i < _tetrahedra.size(); i++)
    {
        const Tetra* tetra = _tetrahedra[i];

        vector<double> coords;
        vector<int> indices = {3, 0, 1, 2, 3, 0, 1, 3, 3, 0, 2, 3, 3, 1, 2, 3};

        for (size_t j = 0; j < 4; j++)
        {
            const Vec* v = tetra->_vertices[j];
            coords.push_back(v->x());
            coords.push_back(v->y());
            coords.push_back(v->z());
        }
        // write the edges of the cell to the plot files
        // if (bounds.zmin() <= 0 && bounds.zmax() >= 0) plotxy.writePolyhedron(coords, indices);
        // if (bounds.ymin() <= 0 && bounds.ymax() >= 0) plotxz.writePolyhedron(coords, indices);
        // if (bounds.xmin() <= 0 && bounds.xmax() >= 0) plotyz.writePolyhedron(coords, indices);
        if (_tetrahedra.size() <= 1000) plotxyz.writePolyhedron(coords, indices);
    }
}

////////////////////////////////////////////////////////////////////

Box TetraMeshSnapshot::extent() const
{
    return _extent;
}

////////////////////////////////////////////////////////////////////

int TetraMeshSnapshot::numEntities() const
{
    return _tetrahedra.size();
}

////////////////////////////////////////////////////////////////////

Position TetraMeshSnapshot::position(int m) const
{
    return _tetrahedra[m]->position();
}

////////////////////////////////////////////////////////////////////

double TetraMeshSnapshot::volume(int m) const
{
    return _tetrahedra[m]->volume();
}

////////////////////////////////////////////////////////////////////

Box TetraMeshSnapshot::extent(int m) const
{
    return _tetrahedra[m]->extent();
}

////////////////////////////////////////////////////////////////////

double TetraMeshSnapshot::density(int m) const
{
    return _rhov[m];
}

////////////////////////////////////////////////////////////////////

double TetraMeshSnapshot::density(Position bfr) const
{
    int m = cellIndex(bfr);
    return m >= 0 ? _rhov[m] : 0;
}

////////////////////////////////////////////////////////////////////

double TetraMeshSnapshot::mass() const
{
    return _mass;
}

////////////////////////////////////////////////////////////////////

Position TetraMeshSnapshot::generatePosition(int m) const
{
    // get loop-invariant information about the cell
    const Box& box = _tetrahedra[m]->extent();

    // generate random points in the enclosing box until one happens to be inside the cell
    for (int i = 0; i < 10000; i++)
    {
        Position r = random()->position(box);
        if (_tetrahedra[m]->inside(r)) return r;
    }
    log()->error("can't find random poisition in tetrahedron");
    return Position();
    // throw FATALERROR("Can't find random position in cell");
}

////////////////////////////////////////////////////////////////////

Position TetraMeshSnapshot::generatePosition() const
{
    // if there are no sites, return the origin
    if (_tetrahedra.empty()) return Position();

    // select a site according to its mass contribution
    int m = NR::locateClip(_cumrhov, random()->uniform());

    return generatePosition(m);
}

////////////////////////////////////////////////////////////////////

int TetraMeshSnapshot::cellIndex(Position bfr) const
{
    int tetras = _tetrahedra.size();
    for (int i = 0; i < tetras; i++)
    {
        const Tetra* tetra = _tetrahedra[i];
        // speed up by rejecting all points that are not in de bounding box
        if (tetra->Box::contains(bfr) && tetra->Tetra::inside(bfr)) return i;
    }
    return -1;
}

////////////////////////////////////////////////////////////////////

const Array& TetraMeshSnapshot::properties(int m) const
{
    return _cells[m]->properties();
}

////////////////////////////////////////////////////////////////////

int TetraMeshSnapshot::nearestEntity(Position bfr) const
{
    return _blocktrees.size() ? cellIndex(bfr) : -1;
}

////////////////////////////////////////////////////////////////////

void TetraMeshSnapshot::getEntities(EntityCollection& entities, Position bfr) const
{
    entities.addSingle(cellIndex(bfr));
}

////////////////////////////////////////////////////////////////////

class TetraMeshSnapshot::MySegmentGenerator : public PathSegmentGenerator
{
    const TetraMeshSnapshot* _grid{nullptr};
    int _mr{-1};
    // true if the path is was ever inside the convex hull
    bool wasInside = false;
    int enteringFace = -1;

public:
    MySegmentGenerator(const TetraMeshSnapshot* grid) : _grid(grid) {}

    std::ofstream out;

    void startwriting()
    {
        if (!out.is_open()) out.open("photon.txt");
    }

    void stopwriting() { out.close(); }

    void write(const Vec& exit, double ds, int face)
    {
        out << "photon=" << _mr << "," << face << "\n";
        out << "ds=" << ds << "\n";
        out << "r=" << r().x() << "," << r().y() << "," << r().z() << "\n";
        out << "exit=" << exit.x() << "," << exit.y() << "," << exit.z() << "\n";
        out << "k=" << k().x() << "," << k().y() << "," << k().z() << "\n";
    }

    bool next() override
    {
        if (state() == State::Unknown)
        {
            startwriting();

            // try moving the photon packet inside the grid; if this is impossible, return an empty path
            if (!moveInside(_grid->extent(), _grid->_eps)) return false;

            // get the index of the cell containing the current position
            _mr = _grid->cellIndex(r());

            if (_mr == -1)
                setState(State::Outside);
            else
                wasInside = true;  // setState has already been called

            // if the photon packet started outside the grid, return the corresponding nonzero-length segment;
            // otherwise fall through to determine the first actual segment
            if (ds() > 0.) return true;

            write(r(), 0, -1);
        }

        // intentionally falls through
        if (state() == State::Inside)
        {
            // loop in case no exit point was found (which should happen only rarely)
            while (true)
            {
                // initialize the smallest nonnegative intersection distance and corresponding index
                double sq = DBL_MAX;       // very large, but not infinity (so that infinite si values are discarded)
                const int NO_INDEX = -99;  // meaningless cell index
                int mq = NO_INDEX;

                // plucker coords
                const Plucker ray = Plucker(r(), k());
                const Tetra* tetra = _grid->_tetrahedra[_mr];

                int leavingFace = -1;
                std::array<double, 3> prods;
                if (enteringFace == -1)  // this should only be done once
                {
                    for (int face = 0; face < 4; face++)
                    {
                        if (tetra->intersects(prods, ray, face, true))
                        {
                            leavingFace = face;
                            break;
                        }
                    }
                }
                else
                {
                    int t0 = (enteringFace + 1) % 4;
                    int t1 = (enteringFace + 2) % 4;
                    int t2 = (enteringFace + 3) % 4;
                    // if face is even we should swap two edges
                    if (enteringFace % 2 == 0) std::swap(t0, t2);

                    int edges[3];
                    edges[0] = tetra->getEdge(t0, enteringFace);
                    edges[1] = tetra->getEdge(t1, enteringFace);
                    edges[2] = tetra->getEdge(t2, enteringFace);

                    // 2 step decision tree
                    int i = 0, prev_i;
                    for (int _ = 0; _ < 2; _++)
                    {
                        int e = edges[i];
                        if (e < 0)
                            prods[i] = -Plucker::dot(ray, *tetra->_edges[-e]);
                        else
                            prods[i] = Plucker::dot(ray, *tetra->_edges[e]);

                        // if clockwise move clockwise
                        prev_i = i;
                        if (prods[i] < 0)
                            i = (i + 1) % 3;
                        else  // if counter move counter
                            i = (i - 1) % 3;
                    }

                    // if 2 clockwise: face=t0
                    // if 2 counter: face=t0
                    // if clockwise then counter: face=t2
                    // if counter then clockwise: face=t1

                    // and calculate last inner products
                    if (prods[0] == prods[i])
                    {
                        leavingFace = t0;

                        
                        
                    }
                    else
                    {
                        if (prods[0] < 0)
                            leavingFace = t2;
                        else
                            leavingFace = t1;
                    }
                }

                Vec exit;
                if (leavingFace != -1)
                {
                    for (int i = 0; i < 3; i++)
                    {
                        exit += tetra->_vertices[leavingFace._indices[i]] * w[i];
                    }
                    sq = (exit - r()).norm();
                    mq = tetra->_neighbors[leavingIndex];
                }

                // if no exit point was found, advance the current point by a small distance,
                // recalculate the cell index, and return to the start of the loop
                if (mq == NO_INDEX)
                {
                    propagater(_grid->_eps);
                    _mr = _grid->cellIndex(r());

                    // if we're outside the domain, terminate the path without returning a path segment
                    if (_mr < 0)
                    {
                        setState(State::Outside);
                        return false;
                    }
                }
                // otherwise set the current point to the exit point and return the path segment
                else
                {
                    propagater(sq + _grid->_eps);
                    setSegment(_mr, sq);
                    write(exit, sq, leavingFace);
                    _mr = mq;

                    // if we're outside the domain, terminate the path after returning this path segment
                    if (_mr < 0) setState(State::Outside);
                    return true;
                }
            }
        }

        if (state() == State::Outside)
        {
            // outside the convex hull and inside extent
            if (wasInside && _grid->extent().contains(r()))
            {
                double t_x, t_y, t_z;
                if (kx() < 0)
                    t_x = (_grid->extent().xmin() - rx()) / kx();
                else
                    t_x = (_grid->extent().xmax() - rx()) / kx();
                if (ky() < 0)
                    t_y = (_grid->extent().ymin() - ry()) / ky();
                else
                    t_y = (_grid->extent().ymax() - ry()) / ky();
                if (kz() < 0)
                    t_z = (_grid->extent().zmin() - rz()) / kz();
                else
                    t_z = (_grid->extent().zmax() - rz()) / kz();

                double ds = min(t_x, min(t_y, t_z));
                propagater(ds + _grid->_eps);
                write(r(), ds, -1);
                stopwriting();
                return false;
            }

            // figure out if the path intersects the convex hull
            if (!wasInside)
            {
                const Plucker ray = Plucker(r(), k());
                int enteringfaces = 0;
                int entering = -1;

                // std::array<double, 3> w;
                // std::array<int, 3> vertices;
                // for (int i = 0; i < _grid->_tetrahedra.size(); i++)
                // {
                //     const Tetra* tetra = _grid->_tetrahedra[i];
                //     for (int f = 0; f < 4; f++)
                //     {
                //         int n = tetra->_neighbors[f];
                //         // edge face
                //         if (n == -1)
                //         {
                //             Face face = Face(tetra, f);
                //             if (face.intersects(w, ray, false))
                //             {
                //                 _mr = i;

                //                 if (enteringfaces++ == 0)
                //                 {
                //                     Vec exit;
                //                     for (int i = 0; i < 3; i++)
                //                     {
                //                         exit += tetra->_vertices[face._indices[i]] * w[i];
                //                     }

                //                     double sq = (exit - r()).norm();
                //                     propagater(sq + _grid->_eps);
                //                     write(exit, sq, f);
                //                 }
                //                 else
                //                 {
                //                     _grid->log()->error("too many entering faces for outside ray");
                //                 }
                //                 // break;
                //             }
                //         }
                //     }
                // }
                if (enteringfaces > 0)
                {
                    wasInside = true;
                    setState(State::Inside);
                    return true;
                }
            }
        }
        stopwriting();
        return false;
    }
};

////////////////////////////////////////////////////////////////////

void TetraMeshSnapshot::getEntities(EntityCollection& entities, Position bfr, Direction bfk) const
{
    // initialize a path segment generator
    MySegmentGenerator generator(this);
    generator.start(bfr, bfk);

    // determine and store the path segments in the entity collection
    entities.clear();
    while (generator.next())
    {
        entities.add(generator.m(), generator.ds());
    }
}

////////////////////////////////////////////////////////////////////

std::unique_ptr<PathSegmentGenerator> TetraMeshSnapshot::createPathSegmentGenerator() const
{
    return std::make_unique<MySegmentGenerator>(this);
}

////////////////////////////////////////////////////////////////////
