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

inline double TetraMeshSnapshot::Plucker::dot(const Plucker& a, const Plucker& b)
{
    return Vec::dot(a.U, b.V) + Vec::dot(b.U, a.V);
}

////////////////////////////////////////////////////////////////////

TetraMeshSnapshot::Edge::Edge(const int i1, const int i2, const Vec& v1, const Vec& v2)
    : Plucker(v1, v2 - v1), i1(i1), i2(i2)
{}

////////////////////////////////////////////////////////////////////

bool TetraMeshSnapshot::Edge::operator==(const Edge& edge) const
{
    return (i1 == edge.i1 && i2 == edge.i2) || (i1 == edge.i2 && i2 == edge.i1);
}

////////////////////////////////////////////////////////////////////

size_t TetraMeshSnapshot::Edge::hash(const int i1, const int i2)
{
    size_t max = std::max(i1, i2);
    size_t min = std::min(i1, i2);
    return (max << 32) + min;
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

    _volume = 1 / 6.
              * abs(Vec::dot(Vec::cross(*_vertices[1] - *_vertices[0], *_vertices[2] - *_vertices[0]),
                             *_vertices[3] - *_vertices[0]));
}

////////////////////////////////////////////////////////////////////

double TetraMeshSnapshot::Tetra::getProd(const Plucker& ray, int t1, int t2) const
{
    double prod = 1;
    if (t1 > t2)
    {
        std::swap(t1, t2);
        prod *= -1;
    }
    int e = (t1 == 0) ? t2 - 1 : t1 + t2;
    Edge* edge = _edges[e];
    if (_vertex_indices[t1] != edge->i1)  // not the same order
    {
        prod *= -1;
    }
    return prod * Plucker::dot(ray, *edge);
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

void TetraMeshSnapshot::Tetra::addEdges(EdgeMap& edgemap)
{

    for (int t1 = 0; t1 < 3; t1++)
    {
        for (int t2 = t1 + 1; t2 < 4; t2++)
        {
            // _edgemap
            int i1 = _vertex_indices[t1];
            int i2 = _vertex_indices[t2];
            size_t key = Edge::hash(i1, i2);
            auto it = edgemap.find(key);

            Edge* edge;
            if (it == edgemap.end())  // not in map
            {
                edge = new Edge(i1, i2, *_vertices[t1], *_vertices[t2]);
                // add this to the map
                edgemap[key] = edge;
            }
            else  // already in map
            {
                edge = it->second;
            }

            // edges are labelled:
            //  0  1  2  3  4  5
            // 01 02 03 12 13 23
            _edges[(t1 == 0) ? t2 - 1 : t1 + t2] = edge;
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

bool TetraMeshSnapshot::Tetra::intersects(std::array<double, 3>& barycoords, const Plucker& ray, int face,
                                          bool leaving) const
{
    std::array<int, 3> t = clockwiseVertices(face);

    double sum = 0;
    for (int i = 0; i < 3; i++)
    {
        // edges: 12, 20, 01
        // verts:  0,  1,  2
        double prod = getProd(ray, t[(i + 1) % 3], t[(i + 2) % 3]);
        if (leaving != (prod < 0)) return false;
        barycoords[i] = prod;
        sum += prod;
    }
    for (int i = 0; i < 3; i++) barycoords[i] /= sum;

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
    // std::array<double, 3> bary = getBary(_vertices);
    // double sum = 0.;
    // for (int i = 0; i < 3; i++)
    // {
    //     if (bary[i] < 0.) return false;
    //     sum += bary[i];
    // }
    // return 1 - sum >= 0;

    // very poor implementation
    const Vec& v0 = *_vertices[0];
    const Vec& v1 = *_vertices[1];
    const Vec& v2 = *_vertices[2];
    const Vec& v3 = *_vertices[3];
    return SameSide(v0, v1, v2, v3, bfr) && SameSide(v1, v2, v3, v0, bfr) && SameSide(v2, v3, v0, v1, bfr)
           && SameSide(v3, v0, v1, v2, bfr);

    // Vec AB = *_vertices[1] - *_vertices[0];
    // Vec AC = *_vertices[2] - *_vertices[0];
    // Vec AD = *_vertices[3] - *_vertices[0];
    // Vec AP = bfr - *_vertices[0];

    // double volume_ABC = Vec::dot(Vec::cross(AB, AC), AD);
    // double volume_PBC = Vec::dot(Vec::cross(AP, AC), AD);
    // double volume_PAC = Vec::dot(Vec::cross(AB, AP), AD);
    // double volume_PAB = Vec::dot(Vec::cross(AB, AC), AP);

    // return (volume_ABC > 0) == (volume_PBC > 0) == (volume_PAC > 0) == (volume_PAB > 0);
}

////////////////////////////////////////////////////////////////////

Vec TetraMeshSnapshot::Tetra::calcExit(const std::array<double, 3>& barycoords, int face) const
{
    std::array<int, 3> t = Tetra::clockwiseVertices(face);
    Vec exit;
    for (int i = 0; i < 3; i++) exit += *_vertices[t[i]] * barycoords[i];
    return exit;
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

std::array<int, 3> TetraMeshSnapshot::Tetra::clockwiseVertices(int face)
{
    std::array<int, 3> t = {(face + 1) % 4, (face + 2) % 4, (face + 3) % 4};
    // if face is even we should swap two edges
    if (face % 2 == 0) std::swap(t[0], t[2]);
    return t;
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

    // std::array<double, 3> getBary(const std::array<Vec*, 4>& vertices)
    // {
    //     double T[3][3];
    //     double T_inv[3][3];

    //     for (int i = 0; i < 3; i++)
    //     {
    //         T[0][i] = vertices[i]->x() - vertices[3]->x();
    //         T[1][i] = vertices[i]->y() - vertices[3]->y();
    //         T[2][i] = vertices[i]->z() - vertices[3]->z();
    //     }

    //     double det = T[0][0] * (T[1][1] * T[2][2] - T[2][1] * T[1][2])
    //                  - T[0][1] * (T[1][0] * T[2][2] - T[1][2] * T[2][0])
    //                  + T[0][2] * (T[1][0] * T[2][1] - T[1][1] * T[2][0]);

    //     double invdet = 1 / det;

    //     T_inv[0][0] = (T[1][1] * T[2][2] - T[2][1] * T[1][2]) * invdet;
    //     T_inv[0][1] = (T[0][2] * T[2][1] - T[0][1] * T[2][2]) * invdet;
    //     T_inv[0][2] = (T[0][1] * T[1][2] - T[0][2] * T[1][1]) * invdet;
    //     T_inv[1][0] = (T[1][2] * T[2][0] - T[1][0] * T[2][2]) * invdet;
    //     T_inv[1][1] = (T[0][0] * T[2][2] - T[0][2] * T[2][0]) * invdet;
    //     T_inv[1][2] = (T[1][0] * T[0][2] - T[0][0] * T[1][2]) * invdet;
    //     T_inv[2][0] = (T[1][0] * T[2][1] - T[2][0] * T[1][1]) * invdet;
    //     T_inv[2][1] = (T[2][0] * T[0][1] - T[0][0] * T[2][1]) * invdet;
    //     T_inv[2][2] = (T[0][0] * T[1][1] - T[1][0] * T[0][1]) * invdet;

    //     std::array<double, 3> bary;
    //     double vec[3] = {vertices[3]->x(), vertices[3]->y(), vertices[3]->z()};
    //     for (int i = 0; i < 3; i++)
    //     {
    //         for (int j = 0; j < 3; j++)
    //         {
    //             bary[i] += T_inv[i][j] * vec[j];
    //         }
    //     }
    //     return bary;
    // }
}

////////////////////////////////////////////////////////////////////

TetraMeshSnapshot::TetraMeshSnapshot() {}

////////////////////////////////////////////////////////////////////

TetraMeshSnapshot::~TetraMeshSnapshot()
{
    for (auto cell : _cells) delete cell;
    for (auto tetra : _tetrahedra) delete tetra;
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
                    tetra->addEdges(_edgemap);

                    _tetrahedra.push_back(tetra);
                }
            }
        }
    }
    _tetrahedra.shrink_to_fit();
    numTetra = _tetrahedra.size();

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

    log()->info("Done computing Delaunay tetrahedralization with " + std::to_string(_tetrahedra.size()) + " cells");
    log()->info("number of edges: " + std::to_string(_edgemap.size()));
    log()->info("bucket count for edgeset: " + std::to_string(_edgemap.bucket_count()));
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

void TetraMeshSnapshot::writeGridPlotFiles(const SimulationItem* /*probe*/) const
{
    std::ofstream outputFile("data/input.txt");
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
    // for (size_t i = 0; i < _tetrahedra.size(); i++)
    // {
    //     outputFile << "circumsphere=" << i << "\n";
    //     outputFile << centers[i].x() << "," << centers[i].y() << "," << centers[i].z() << "\n";
    //     outputFile << radii[i] << "\n";
    // }

    outputFile.close();
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
    // https://vcg.isti.cnr.it/activities/OLD/geometryegraphics/pointintetraedro.html
    double s = random()->uniform();
    double t = random()->uniform();
    double u = random()->uniform();
    if (s + t > 1.0)
    {  // cut'n fold the cube into a prism

        s = 1.0 - s;
        t = 1.0 - t;
    }
    if (t + u > 1.0)
    {  // cut'n fold the prism into a tetrahedron

        double tmp = u;
        u = 1.0 - s - t;
        t = 1.0 - tmp;
    }
    else if (s + t + u > 1.0)
    {

        double tmp = u;
        u = s + t + u - 1.0;
        s = 1 - t - tmp;
    }
    double a = 1 - s - t - u;  // a,s,t,u are the barycentric coordinates of the random point.

    const auto& vert = _tetrahedra[m]->_vertices;
    return Position(*vert[0] * a + *vert[1] * s + *vert[2] * t + *vert[3] * u);
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
    for (int i = 0; i < numTetra; i++)
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
    return cellIndex(bfr);
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
    bool wasInside = false;  // true if the path is was ever inside the convex hull
    int enteringFace = -1;

public:
    MySegmentGenerator(const TetraMeshSnapshot* grid) : _grid(grid) {}

    bool next() override
    {
        if (state() == State::Unknown)
        {
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
        }

        // intentionally falls through
        if (state() == State::Inside)
        {
            // loop in case no exit point was found (which should happen only rarely)
            while (true)
            {
                const Plucker ray = Plucker(r(), k());
                const Tetra* tetra = _grid->_tetrahedra[_mr];

                // temp variable to store plucker products and eventually the barycentric coordinates
                std::array<double, 3> prods;

                int leavingFace = -1;
                if (enteringFace == -1)  // this should only be done once max per ray
                {
                    for (int face = 0; face < 4; face++)
                    {
                        if (tetra->intersects(prods, ray, face, true))
                        {
                            leavingFace = face;
                            break;
                        }
                    }
                    // if (leavingFace == -1)
                    // {
                    //     _grid->log()->warning("no leaving face found!");
                    // }
                }
                else
                {
                    std::array<int, 3> t = Tetra::clockwiseVertices(enteringFace);

                    // 2 step decision tree
                    prods[0] = tetra->getProd(ray, t[0], enteringFace);
                    bool clockwise0 = prods[0] < 0;
                    int i = clockwise0 ? 1 : 2;  // if (counter)clockwise move (counter)clockwise
                    prods[i] = tetra->getProd(ray, t[i], enteringFace);

                    // if 2 clockwise: face=t0
                    // if 2 counter: face=t0
                    // if clockwise then counter: face=t2
                    // if counter then clockwise: face=t1

                    if (clockwise0 == (prods[i] < 0))
                        leavingFace = t[0];
                    else if (clockwise0)
                        leavingFace = t[2];
                    else
                        leavingFace = t[1];

                    // get prods (this calculates one prod too many but is much cleaner)
                    tetra->intersects(prods, ray, leavingFace, true);
                }

                // if no exit point was found, advance the current point by a small distance,
                // recalculate the cell index, and return to the start of the loop
                if (leavingFace == -1)
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
                    Vec exit = tetra->calcExit(prods, leavingFace);

                    double ds = (exit - r()).norm();
                    int next_mr = tetra->_neighbors[leavingFace];
                    // we could assert if r+exit and r+sq*k are the same

                    propagater(ds + _grid->_eps);
                    setSegment(_mr, ds);

                    // set enteringFace if there is a neighboring cell
                    if (next_mr != -1)
                    {
                        auto& neighbors = _grid->_tetrahedra[next_mr]->_neighbors;
                        enteringFace =
                            std::distance(neighbors.begin(), std::find(neighbors.begin(), neighbors.end(), _mr));
                    }
                    else
                    {
                        enteringFace = -1;
                    }
                    // set new cell
                    _mr = next_mr;

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

                double ds = min({t_x, t_y, t_z});
                propagater(ds + _grid->_eps);
                setSegment(_mr, ds);
                return false;
            }

            // figure out if the path intersects the convex hull
            if (!wasInside)
            {
                const Plucker ray = Plucker(r(), k());
                std::array<double, 3> barycoords;
                for (int i = 0; i < _grid->numTetra; i++)
                {
                    const Tetra* tetra = _grid->_tetrahedra[i];
                    for (int face = 0; face < 4; face++)
                    {
                        // convex hull face
                        if (tetra->_neighbors[face] == -1)
                        {
                            if (tetra->intersects(barycoords, ray, face, false))
                            {
                                _mr = i;
                                enteringFace = face;

                                Vec exit = tetra->calcExit(barycoords, enteringFace);

                                double ds = (exit - r()).norm();
                                propagater(ds + _grid->_eps);
                                // setSegment(-1, ds);  // not sure if this is needed

                                wasInside = true;
                                setState(State::Inside);
                                return true;
                            }
                        }
                    }
                }
            }
        }
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
