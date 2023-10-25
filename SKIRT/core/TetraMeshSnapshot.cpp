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
#include <unordered_set>
#include "container.hh"

////////////////////////////////////////////////////////////////////

namespace
{
    // classes used for serializing/deserializing Tetra cell geometry when communicating the results of
    // grid construction between multiple processes with the ProcessManager::broadcastAllToAll() function

    // decorates a std::vector with functions to write serialized versions of various data types
    class SerializedWrite
    {
    private:
        vector<double>& _data;

    public:
        SerializedWrite(vector<double>& data) : _data(data) { _data.clear(); }
        void write(double v) { _data.push_back(v); }
        void write(Vec v) { _data.insert(_data.end(), {v.x(), v.y(), v.z()}); }
        void write(Box v) { _data.insert(_data.end(), {v.xmin(), v.ymin(), v.zmin(), v.xmax(), v.ymax(), v.zmax()}); }
        void write(const vector<int>& v)
        {
            _data.push_back(v.size());
            _data.insert(_data.end(), v.begin(), v.end());
        }
    };

    // decorates a std::vector with functions to read serialized versions of various data types
    class SerializedRead
    {
    private:
        const double* _data;
        const double* _end;

    public:
        SerializedRead(const vector<double>& data) : _data(data.data()), _end(data.data() + data.size()) {}
        bool empty() { return _data == _end; }
        int readInt() { return *_data++; }
        void read(double& v) { v = *_data++; }
        void read(Vec& v)
        {
            v.set(*_data, *(_data + 1), *(_data + 2));
            _data += 3;
        }
        void read(Box& v)
        {
            v = Box(*_data, *(_data + 1), *(_data + 2), *(_data + 3), *(_data + 4), *(_data + 5));
            _data += 6;
        }
        void read(vector<int>& v)
        {
            int n = *_data++;
            v.clear();
            v.reserve(n);
            for (int i = 0; i != n; ++i) v.push_back(*_data++);
        }
    };
}

////////////////////////////////////////////////////////////////////

// class to hold the information about a Tetra cell that is relevant for calculating paths and densities
class TetraMeshSnapshot::Cell : public Box  // enclosing box
{
public:
    Vec _r;                  // site position
    Vec _c;                  // centroid position
    double _volume{0.};      // volume
    vector<int> _neighbors;  // list of neighbor indices in _cells vector
    Array _properties;       // user-defined properties, if any

public:
    // constructor stores the specified site position; the other data members are set to zero or empty
    Cell(Vec r) : _r(r) {}

    // constructor derives the site position from the first three property values and stores the user properties;
    // the other data members are set to zero or empty
    Cell(const Array& prop) : _r(prop[0], prop[1], prop[2]), _properties{prop} {}

    // adjusts the site position with the specified offset
    void relax(double cx, double cy, double cz) { _r += Vec(cx, cy, cz); }

    // initializes the receiver with information taken from the specified fully computed Tetra cell
    void init(voro::voronoicell_neighbor& cell)
    {
        // copy basic geometric info
        double cx, cy, cz;
        cell.centroid(cx, cy, cz);
        _c = Vec(cx, cy, cz) + _r;
        _volume = cell.volume();

        // get the minimal and maximal coordinates of the box enclosing the cell
        vector<double> coords;
        cell.vertices(_r.x(), _r.y(), _r.z(), coords);
        double xmin = DBL_MAX;
        double ymin = DBL_MAX;
        double zmin = DBL_MAX;
        double xmax = -DBL_MAX;
        double ymax = -DBL_MAX;
        double zmax = -DBL_MAX;
        int n = coords.size();
        for (int i = 0; i < n; i += 3)
        {
            xmin = min(xmin, coords[i]);
            ymin = min(ymin, coords[i + 1]);
            zmin = min(zmin, coords[i + 2]);
            xmax = max(xmax, coords[i]);
            ymax = max(ymax, coords[i + 1]);
            zmax = max(zmax, coords[i + 2]);
        }

        // set our inherited Box to this bounding box
        setExtent(xmin, ymin, zmin, xmax, ymax, zmax);

        // copy a list of neighboring cell/site ids
        cell.neighbors(_neighbors);
    }

    // clears the information taken from a Tetra cell so it can be reinitialized
    void clear()
    {
        _volume = 0.;
        _neighbors.clear();
    }

    // initializes the receiver with the volume calculated from imported information
    void init(double volume) { _volume = volume; }

    // returns the cell's site position
    Vec position() const { return _r; }

    // returns the x coordinate of the cell's site position
    double x() const { return _r.x(); }

    // returns the squared distance from the cell's site to the specified point
    double squaredDistanceTo(Vec r) const { return (r - _r).norm2(); }

    // returns the central position in the cell
    Vec centroid() const { return _c; }

    // returns the volume of the cell; overriding volume() function of Box bas class
    double volume() const { return _volume; }

    // returns a list of neighboring cell/site ids
    const vector<int>& neighbors() { return _neighbors; }

    // returns the cell/site user properties, if any
    const Array& properties() { return _properties; }

    // writes the Tetra cell geometry to the serialized data buffer, preceded by the specified cell index,
    // if the cell geometry has been calculated for this cell; otherwise does nothing
    void writeGeometryIfPresent(SerializedWrite& wdata, int m)
    {
        if (!_neighbors.empty())
        {
            wdata.write(m);
            wdata.write(extent());
            wdata.write(_c);
            wdata.write(_volume);
            wdata.write(_neighbors);
        }
    }

    // reads the Tetra cell geometry from the serialized data buffer
    void readGeometry(SerializedRead& rdata)
    {
        rdata.read(*this);  // extent
        rdata.read(_c);
        rdata.read(_volume);
        rdata.read(_neighbors);
    }
};

class TetraMeshSnapshot::Tetra : public Box
{
public:
    Tetra(const vector<Cell*>& _cells, int i, int j, int k, int l)
    {
        _vertices[0] = _cells[i]->_r;
        _vertices[1] = _cells[j]->_r;
        _vertices[2] = _cells[k]->_r;
        _vertices[3] = _cells[l]->_r;

        double xmin = DBL_MAX;
        double ymin = DBL_MAX;
        double zmin = DBL_MAX;
        double xmax = -DBL_MAX;
        double ymax = -DBL_MAX;
        double zmax = -DBL_MAX;
        for (int i = 0; i < 4; i += 3)
        {
            xmin = min(xmin, _vertices[i].x());
            ymin = min(ymin, _vertices[i].y());
            zmin = min(zmin, _vertices[i].z());
            xmax = max(xmax, _vertices[i].x());
            ymax = max(ymax, _vertices[i].y());
            zmax = max(zmax, _vertices[i].z());
        }
        setExtent(Box(xmin, ymin, zmin, xmax, ymax, zmax));

        _vertex_indices[0] = i;
        _vertex_indices[1] = j;
        _vertex_indices[2] = k;
        _vertex_indices[3] = l;
    }

    // compute normal for face 3 (in clockwise direction of vertices 012) and dot with vertex 3
    double orient()
    {
        const Vec e01 = _vertices[1] - _vertices[0];
        const Vec e02 = _vertices[2] - _vertices[0];
        const Vec e03 = _vertices[3] - _vertices[0];
        // this convention makes edges go clockwise around leaving rays from inside the tetrahedron
        double orientation = Vec::dot(Vec::cross(e02, e01), e03);
        if (orientation < 0)
        {
            std::swap(_vertices[0], _vertices[1]);
            std::swap(_vertex_indices[0], _vertex_indices[1]);
        }
        return orientation;
    }

    // return -1 if no shared face
    // return 0-3 for opposite vertex of shared face
    int shareFace(const Tetra* other) const
    {
        int equal_vert = 0;
        int opposite = 0 + 1 + 2 + 3;

        for (size_t i = 0; i < 4; i++)
        {
            const int vi = _vertex_indices[i];

            for (size_t j = 0; j < 4; j++)
            {
                const int vj = other->_vertex_indices[i];

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

    bool equals(const Tetra* other) const
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

    bool SameSide(const Vec& v0, const Vec& v1, const Vec& v2, const Vec& v3, const Vec& pos) const
    {
        Vec normal = Vec::cross(v1 - v0, v2 - v0);
        double dotV4 = Vec::dot(normal, v3 - v0);
        double dotP = Vec::dot(normal, pos - v0);
        return (dotV4 > 0) == (dotP > 0);
    }

    bool inside(const Position& bfr) const
    {
        // very poor implementation
        const Vec& v0 = _vertices[0];
        const Vec& v1 = _vertices[1];
        const Vec& v2 = _vertices[2];
        const Vec& v3 = _vertices[3];
        return SameSide(v0, v1, v2, v3, bfr) && SameSide(v1, v2, v3, v0, bfr) && SameSide(v2, v3, v0, v1, bfr)
               && SameSide(v3, v0, v1, v2, bfr);
    }

    std::array<Vec, 4> _vertices;
    std::array<int, 4> _neighbors = {-1, -1, -1, -1};
    std::array<int, 4> _vertex_indices;
};

class TetraMeshSnapshot::Plucker
{
    Vec U, V;

public:
    Plucker() {}

    Plucker(const Vec& U, const Vec& V) : U(U), V(V) {}

    static inline Plucker createFromDir(const Vec& pos, const Vec& dir) { return Plucker(dir, Vec::cross(dir, pos)); }

    static inline Plucker createFromVertices(const Vec& v1, const Vec& v2)
    {
        Plucker p;
        p.U = v2 - v1;
        p.V = Vec::cross(p.U, v1);
        return p;
    }

    // give plucker coord of all edges of a certain face
    static inline void face(std::array<Plucker, 3>& edges, const std::array<Vec, 4>& vertices, int face)
    {
        int v1 = (face + 1) % 4;
        int v2 = (face + 2) % 4;
        int v3 = (face + 3) % 4;

        // if face is even we should swap two edges
        if (face % 2 == 0) std::swap(v2, v3);

        edges[0] = Plucker(vertices[v1], vertices[v2]);
        edges[1] = Plucker(vertices[v2], vertices[v3]);
        edges[2] = Plucker(vertices[v3], vertices[v1]);
    }

    // permuted inner product
    static inline double dot(const Plucker& a, const Plucker& b) { return Vec::dot(a.U, b.V) + Vec::dot(b.U, a.V); }
};

////////////////////////////////////////////////////////////////////

namespace
{
    // function to compare two points according to the specified axis (0,1,2)
    bool lessthan(Vec p1, Vec p2, int axis)
    {
        switch (axis)
        {
            case 0:  // split on x
                if (p1.x() < p2.x()) return true;
                if (p1.x() > p2.x()) return false;
                if (p1.y() < p2.y()) return true;
                if (p1.y() > p2.y()) return false;
                if (p1.z() < p2.z()) return true;
                return false;
            case 1:  // split on y
                if (p1.y() < p2.y()) return true;
                if (p1.y() > p2.y()) return false;
                if (p1.z() < p2.z()) return true;
                if (p1.z() > p2.z()) return false;
                if (p1.x() < p2.x()) return true;
                return false;
            case 2:  // split on z
                if (p1.z() < p2.z()) return true;
                if (p1.z() > p2.z()) return false;
                if (p1.x() < p2.x()) return true;
                if (p1.x() > p2.x()) return false;
                if (p1.y() < p2.y()) return true;
                return false;
            default:  // this should never happen
                return false;
        }
    }

    template<typename T> bool invec(const vector<T>& vec, const T& e)
    {
        return std::find(vec.begin(), vec.end(), e) != vec.end();
    }

    double det4(double x0, double y0, double z0, double x1, double y1, double z1, double x2, double y2, double z2,
                double x3, double y3, double z3)
    {
        // return x0 * (y1 * (z2 - z3) - y2 * (z1 - z3) + y3 * (z1 - z2))
        //        - y0 * (x1 * (z2 - z3) - x2 * (z1 - z3) + x3 * (z1 - z2))
        //        + z0 * (x1 * (y2 - y3) - x2 * (y1 - y3) + x3 * (y1 - y2));

        return x0 * y1 * z2 - x0 * y1 * z3 - x0 * y2 * z1 + x0 * y2 * z3 + x0 * y3 * z1 - x0 * y3 * z2 - x1 * y0 * z2
               + x1 * y0 * z3 + x1 * y2 * z0 - x1 * y2 * z3 - x1 * y3 * z0 + x1 * y3 * z2 + x2 * y0 * z1 - x2 * y0 * z3
               - x2 * y1 * z0 + x2 * y1 * z3 + x2 * y3 * z0 - x2 * y3 * z1 - x3 * y0 * z1 + x3 * y0 * z2 + x3 * y1 * z0
               - x3 * y1 * z2 - x3 * y2 * z0 + x3 * y2 * z1;
    }
}

////////////////////////////////////////////////////////////////////

// class to hold a node in the binary search tree (see en.wikipedia.org/wiki/Kd-tree)
class TetraMeshSnapshot::Node
{
private:
    int _m;        // index in _cells to the site defining the split at this node
    int _axis;     // split axis for this node (0,1,2)
    Node* _up;     // ptr to the parent node
    Node* _left;   // ptr to the left child node
    Node* _right;  // ptr to the right child node

    // returns the square of its argument
    static double sqr(double x) { return x * x; }

public:
    // constructor stores the specified site index and child pointers (which may be null)
    Node(int m, int depth, Node* left, Node* right) : _m(m), _axis(depth % 3), _up(0), _left(left), _right(right)
    {
        if (_left) _left->setParent(this);
        if (_right) _right->setParent(this);
    }

    // destructor destroys the children
    ~Node()
    {
        delete _left;
        delete _right;
    }

    // sets parent pointer (called from parent's constructor)
    void setParent(Node* up) { _up = up; }

    // returns the corresponding data member
    int m() const { return _m; }
    Node* up() const { return _up; }
    Node* left() const { return _left; }
    Node* right() const { return _right; }

    // returns the apropriate child for the specified query point
    Node* child(Vec bfr, const vector<Cell*>& cells) const
    {
        return lessthan(bfr, cells[_m]->position(), _axis) ? _left : _right;
    }

    // returns the other child than the one that would be apropriate for the specified query point
    Node* otherChild(Vec bfr, const vector<Cell*>& cells) const
    {
        return lessthan(bfr, cells[_m]->position(), _axis) ? _right : _left;
    }

    // returns the squared distance from the query point to the split plane
    double squaredDistanceToSplitPlane(Vec bfr, const vector<Cell*>& cells) const
    {
        switch (_axis)
        {
            case 0:  // split on x
                return sqr(cells[_m]->position().x() - bfr.x());
            case 1:  // split on y
                return sqr(cells[_m]->position().y() - bfr.y());
            case 2:  // split on z
                return sqr(cells[_m]->position().z() - bfr.z());
            default:  // this should never happen
                return 0;
        }
    }

    // returns the node in this subtree that represents the site nearest to the query point
    Node* nearest(Vec bfr, const vector<Cell*>& cells)
    {
        // recursively descend the tree until a leaf node is reached, going left or right depending on
        // whether the specified point is less than or greater than the current node in the split dimension
        Node* current = this;
        while (Node* child = current->child(bfr, cells)) current = child;

        // unwind the recursion, looking for the nearest node while climbing up
        Node* best = current;
        double bestSD = cells[best->m()]->squaredDistanceTo(bfr);
        while (true)
        {
            // if the current node is closer than the current best, then it becomes the current best
            double currentSD = cells[current->m()]->squaredDistanceTo(bfr);
            if (currentSD < bestSD)
            {
                best = current;
                bestSD = currentSD;
            }

            // if there could be points on the other side of the splitting plane for the current node
            // that are closer to the search point than the current best, then ...
            double splitSD = current->squaredDistanceToSplitPlane(bfr, cells);
            if (splitSD < bestSD)
            {
                // move down the other branch of the tree from the current node looking for closer points,
                // following the same recursive process as the entire search
                Node* other = current->otherChild(bfr, cells);
                if (other)
                {
                    Node* otherBest = other->nearest(bfr, cells);
                    double otherBestSD = cells[otherBest->m()]->squaredDistanceTo(bfr);
                    if (otherBestSD < bestSD)
                    {
                        best = otherBest;
                        bestSD = otherBestSD;
                    }
                }
            }

            // move up to the parent until we meet the top node
            if (current == this) break;
            current = current->up();
        }
        return best;
    }
};

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

void TetraMeshSnapshot::readAndClose()
{
    // read the site info into memory
    Array prop;
    while (infile()->readRow(prop)) _cells.push_back(new Cell(prop));

    // close the file
    Snapshot::readAndClose();

    // if we are allowed to build a Tetra mesh
    if (!_foregoTetraMesh)
    {
        // calculate the Tetra cells
        buildMesh(false);

        // if a mass density policy has been set, calculate masses and densities and build the search data structure
        if (hasMassDensityPolicy()) calculateDensityAndMass();
        if (hasMassDensityPolicy() || needGetEntities()) buildSearchPerBlock();
    }

    // if we forego building a Tetra mesh, there is a density policy by definition
    else
    {
        calculateVolume();
        calculateDensityAndMass();
        buildSearchSingle();
    }
}

////////////////////////////////////////////////////////////////////

void TetraMeshSnapshot::setExtent(const Box& extent)
{
    _extent = extent;
    _eps = 1e-12 * extent.widths().norm();
}

////////////////////////////////////////////////////////////////////

void TetraMeshSnapshot::foregoTetraMesh()
{
    _foregoTetraMesh = true;
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
    buildSearchPerBlock();
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
    buildSearchPerBlock();
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
    buildSearchPerBlock();
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
    ///////////////////////////////////////////////////////////////////////////
    setExtent(Box(-1, -1, -1, 1, 1, 1));

    vector<Vec> sites;
    // for (size_t i = 0; sites.size() < 15; i++)
    // {
    //     Position r = random()->position(_extent);
    //     sites.push_back(r);
    // }

    // sites.push_back(Vec( 0.2,  0.2,  0));
    // sites.push_back(Vec( 0.2, -0.2,  0));
    sites.push_back(Vec(0.4, 0.4, 0));
    sites.push_back(Vec(0.4, -0.4, 0));

    sites.push_back(Vec(-0.5, 0, 0));
    sites.push_back(Vec(0, 0, 0.5));
    sites.push_back(Vec(0, 0, -0.5));

    // prepare the data
    int n = sites.size();
    _cells.resize(n);
    for (int m = 0; m != n; ++m) _cells[m] = new Cell(sites[m]);

    const double A = 100.;

    ///////////////////////////////////////////////////////////////////////////

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

        // add the retained original sites to a temporary Tetra container, using the cell index m as ID
        voro::container vcon(-A, A, -A, A, -A, A, _nb, _nb, _nb, false, false, false, 16);
        for (int m = 0; m != numCells; ++m)
        {
            Vec r = _cells[m]->position();
            vcon.put(m, r.x(), r.y(), r.z());
        }

        // compute the cell in the Tetra tesselation corresponding to each site
        // and store the cell's centroid (relative to the site position) as the relaxation offset
        log()->info("Relaxing Tetra tessellation with " + std::to_string(numCells) + " cells");
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
                        if (numDone == 0) log()->infoIfElapsed("Computed Tetra cells: ", logProgressChunkSize);
                    }
                } while (vloop.inc());
            if (numDone > 0) log()->infoIfElapsed("Computed Tetra cells: ", numDone);
        });

        // communicate the calculated offsets between parallel processes, if needed, and apply them to the cells
        ProcessManager::sumToAll(offsets.data());
        for (int m = 0; m != numCells; ++m) _cells[m]->relax(offsets(m, 0), offsets(m, 1), offsets(m, 2));
    }

    // ========= FINAL GRID =========

    // repeat grid construction until none of the cells have zero volume
    int numIterations = 0;
    while (true)
    {
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
        log()->info("Constructing Tetra tessellation with " + std::to_string(numCells) + " cells");
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

        ////////////////////////////////////////////////////////////////////////////////////////////
        for (int i = 0; i < numCells; i++)
        {
            Cell* c1 = _cells[i];

            for (int j : c1->neighbors())
            {
                if (j < 0) continue;

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

                        double r0 = x0 * x0 + y0 * y0 + z0 * z0;
                        double r1 = x1 * x1 + y1 * y1 + z1 * z1;
                        double r2 = x2 * x2 + y2 * y2 + z2 * z2;
                        double r3 = x3 * x3 + y3 * y3 + z3 * z3;

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
                                    break;
                                }
                            }
                            if (!delaunay) break;
                        }
                        if (!delaunay) continue;

                        // orient tetrahedron in the same consistent way
                        tetra->orient();

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

        ////////////////////////////////////////////////////////////////////////////////////////////

        // discover invalid cells with zero volume and/or with neighbors that are not mutual
        log()->info("Verifying Tetra tessellation");
        std::set<int> invalid;
        for (int m = 0; m < numCells; m++)
        {
            if (!_cells[m]->volume()) invalid.insert(m);
            for (int m1 : _cells[m]->neighbors())
            {
                if (m1 >= 0)
                {
                    const vector<int>& neighbors1 = _cells[m1]->neighbors();
                    if (std::find(neighbors1.begin(), neighbors1.end(), m) == neighbors1.end())
                    {
                        invalid.insert(m);
                        invalid.insert(m1);
                    }
                }
            }
        }

        // break from loop if no invalid cells were found
        if (invalid.empty()) break;

        // give up after a given number of iterations
        if (++numIterations == maxConstructionIterations)
        {
            throw FATALERROR("Still " + std::to_string(invalid.size()) + " invalid Tetra cells after "
                             + std::to_string(maxConstructionIterations)
                             + " iterations of constructing the tessellation");
        }

        // remove invalid cells and prepare to repeat construction
        log()->warning("Removing sites for " + std::to_string(invalid.size())
                       + " invalid Tetra cells and reconstructing the tessellation");
        for (auto it = invalid.begin(); it != invalid.end(); ++it)
        {
            int m = *it;
            delete _cells[m];
            _cells[m] = 0;
        }
        numCells = eraseNullPointers(_cells);
        for (int m = 0; m != numCells; ++m) _cells[m]->clear();
    }

    // ========= STATISTICS =========

    // compile neighbor statistics
    int minNeighbors = INT_MAX;
    int maxNeighbors = 0;
    int64_t totNeighbors = 0;
    for (int m = 0; m < numCells; m++)
    {
        int ns = _cells[m]->neighbors().size();
        totNeighbors += ns;
        minNeighbors = min(minNeighbors, ns);
        maxNeighbors = max(maxNeighbors, ns);
    }
    double avgNeighbors = double(totNeighbors) / numCells;

    // log neighbor statistics
    log()->info("Done computing Tetra tessellation with " + std::to_string(numCells) + " cells");
    log()->info("  Average number of neighbors per cell: " + StringUtils::toString(avgNeighbors, 'f', 1));
    log()->info("  Minimum number of neighbors per cell: " + std::to_string(minNeighbors));
    log()->info("  Maximum number of neighbors per cell: " + std::to_string(maxNeighbors));
}

////////////////////////////////////////////////////////////////////

void TetraMeshSnapshot::calculateVolume()
{
    int numCells = _cells.size();
    for (int m = 0; m != numCells; ++m)
    {
        const Array& prop = _cells[m]->properties();
        double volume = prop[densityIndex()] > 0. ? prop[massIndex()] / prop[densityIndex()] : 0.;
        _cells[m]->init(volume);
    }
}

////////////////////////////////////////////////////////////////////

void TetraMeshSnapshot::calculateDensityAndMass()
{
    // allocate vectors for mass and density
    int numCells = _cells.size();
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
        const Array& prop = _cells[m]->properties();

        // original mass is zero if temperature is above cutoff or if imported mass/density is not positive
        double originalDensity = 0.;
        double originalMass = 0.;
        if (maxT && prop[temperatureIndex()] > maxT)
        {
            numIgnored++;
        }
        else
        {
            double volume = _cells[m]->volume();
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

TetraMeshSnapshot::Node* TetraMeshSnapshot::buildTree(vector<int>::iterator first, vector<int>::iterator last,
                                                      int depth) const
{
    auto length = last - first;
    if (length > 0)
    {
        auto median = length >> 1;
        std::nth_element(first, first + median, last, [this, depth](int m1, int m2) {
            return m1 != m2 && lessthan(_cells[m1]->position(), _cells[m2]->position(), depth % 3);
        });
        return new TetraMeshSnapshot::Node(*(first + median), depth, buildTree(first, first + median, depth + 1),
                                           buildTree(first + median + 1, last, depth + 1));
    }
    return nullptr;
}

////////////////////////////////////////////////////////////////////

void TetraMeshSnapshot::buildSearchPerBlock()
{
    // abort if there are no cells
    int numCells = _cells.size();
    if (!numCells) return;

    log()->info("Building data structures to accelerate searching the Tetra tesselation");

    // -------------  block lists  -------------

    // initialize a vector of nb x nb x nb lists, each containing the cells overlapping a certain block in the domain
    _blocklists.resize(_nb3);

    // add the cell object to the lists for all blocks it may overlap
    int i1, j1, k1, i2, j2, k2;
    for (int m = 0; m != numCells; ++m)
    {
        _extent.cellIndices(i1, j1, k1, _cells[m]->rmin() - Vec(_eps, _eps, _eps), _nb, _nb, _nb);
        _extent.cellIndices(i2, j2, k2, _cells[m]->rmax() + Vec(_eps, _eps, _eps), _nb, _nb, _nb);
        for (int i = i1; i <= i2; i++)
            for (int j = j1; j <= j2; j++)
                for (int k = k1; k <= k2; k++) _blocklists[i * _nb2 + j * _nb + k].push_back(m);
    }

    // compile block list statistics
    int minRefsPerBlock = INT_MAX;
    int maxRefsPerBlock = 0;
    int64_t totalBlockRefs = 0;
    for (int b = 0; b < _nb3; b++)
    {
        int refs = _blocklists[b].size();
        totalBlockRefs += refs;
        minRefsPerBlock = min(minRefsPerBlock, refs);
        maxRefsPerBlock = max(maxRefsPerBlock, refs);
    }
    double avgRefsPerBlock = double(totalBlockRefs) / _nb3;

    // log block list statistics
    log()->info("  Number of blocks in search grid: " + std::to_string(_nb3) + " (" + std::to_string(_nb) + "^3)");
    log()->info("  Average number of cells per block: " + StringUtils::toString(avgRefsPerBlock, 'f', 1));
    log()->info("  Minimum number of cells per block: " + std::to_string(minRefsPerBlock));
    log()->info("  Maximum number of cells per block: " + std::to_string(maxRefsPerBlock));

    // -------------  search trees  -------------

    // for each block that contains more than a predefined number of cells,
    // construct a search tree on the site locations of the cells
    _blocktrees.resize(_nb3);
    for (int b = 0; b < _nb3; b++)
    {
        vector<int>& ids = _blocklists[b];
        if (ids.size() > 9) _blocktrees[b] = buildTree(ids.begin(), ids.end(), 0);
    }

    // compile and log search tree statistics
    int numTrees = 0;
    for (int b = 0; b < _nb3; b++)
        if (_blocktrees[b]) numTrees++;
    log()->info("  Number of search trees: " + std::to_string(numTrees) + " ("
                + StringUtils::toString(100. * numTrees / _nb3, 'f', 1) + "% of blocks)");
}

////////////////////////////////////////////////////////////////////

void TetraMeshSnapshot::buildSearchSingle()
{
    // log the number of sites
    int numCells = _cells.size();
    log()->info("  Number of sites: " + std::to_string(numCells));

    // abort if there are no cells
    if (!numCells) return;

    // construct a single search tree on the site locations of all cells
    log()->info("Building data structure to accelerate searching " + std::to_string(numCells) + " Tetra sites");
    _blocktrees.resize(1);
    vector<int> ids(numCells);
    for (int m = 0; m != numCells; ++m) ids[m] = m;
    _blocktrees[0] = buildTree(ids.begin(), ids.end(), 0);
}

////////////////////////////////////////////////////////////////////

bool TetraMeshSnapshot::isPointClosestTo(Vec r, int m, const vector<int>& ids) const
{
    double target = _cells[m]->squaredDistanceTo(r);
    for (int id : ids)
    {
        if (id >= 0 && _cells[id]->squaredDistanceTo(r) < target) return false;
    }
    return true;
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
            if (n >= 0)
            {
                outputFile << " " << n << ",";
            }
        }
        outputFile << "\n";
    }

    outputFile << "tetrahedra=" << _tetrahedra.size() << "\n";
    for (size_t i = 0; i < _tetrahedra.size(); i++)
    {
        Tetra* tetra = _tetrahedra[i];
        outputFile << "tetrahedron=" << i << "\nvertices=";
        for (size_t l = 0; l < 4; l++)
        {
            outputFile << " " << tetra->_vertex_indices[l] << ",";
        }

        outputFile << "\n" << tetra->_neighbors.size() << " neighbors=";
        for (size_t j = 0; j < 4; j++)
        {
            outputFile << " " << tetra->_neighbors[j] << ",";

            double x = std::round(tetra->_vertices[j].x());
            double y = std::round(tetra->_vertices[j].y());
            double z = std::round(tetra->_vertices[j].z());
            tetra->_vertices[j].set(x, y, z);
        }
        outputFile << "\n";
    }
    for (size_t i = 0; i < _tetrahedra.size(); i++)
    {
        Tetra* tetra = _tetrahedra[i];
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
            const Vec& v = tetra->_vertices[j];
            coords.push_back(v.x());
            coords.push_back(v.y());
            coords.push_back(v.z());
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
    return _cells.size();
}

////////////////////////////////////////////////////////////////////

Position TetraMeshSnapshot::position(int m) const
{
    return Position(_cells[m]->position());
}

////////////////////////////////////////////////////////////////////

Position TetraMeshSnapshot::centroidPosition(int m) const
{
    return Position(_cells[m]->centroid());
}

////////////////////////////////////////////////////////////////////

double TetraMeshSnapshot::volume(int m) const
{
    return _cells[m]->volume();
}

////////////////////////////////////////////////////////////////////

Box TetraMeshSnapshot::extent(int m) const
{
    return _cells[m]->extent();
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
    const Box& box = _cells[m]->extent();
    const vector<int>& neighbors = _cells[m]->neighbors();

    // generate random points in the enclosing box until one happens to be inside the cell
    for (int i = 0; i < 10000; i++)
    {
        Position r = random()->position(box);
        if (isPointClosestTo(r, m, neighbors)) return r;
    }
    throw FATALERROR("Can't find random position in cell");
}

////////////////////////////////////////////////////////////////////

Position TetraMeshSnapshot::generatePosition() const
{
    // if there are no sites, return the origin
    if (_cells.empty()) return Position();

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
    // add entering face index for next

public:
    MySegmentGenerator(const TetraMeshSnapshot* grid) : _grid(grid) {}

    bool next() override
    {
        switch (state())
        {
            case State::Unknown:
            {
                // try moving the photon packet inside the grid; if this is impossible, return an empty path
                if (!moveInside(_grid->extent(), _grid->_eps)) return false;

                // get the index of the cell containing the current position
                _mr = _grid->cellIndex(r());

                // if the photon packet started outside the grid, return the corresponding nonzero-length segment;
                // otherwise fall through to determine the first actual segment
                if (ds() > 0.) return true;
            }

            // intentionally falls through
            case State::Inside:
            {
                // loop in case no exit point was found (which should happen only rarely)
                while (true)
                {
                    // get the site position for this cell
                    Vec pr = _grid->_cells[_mr]->position();

                    // initialize the smallest nonnegative intersection distance and corresponding index
                    double sq = DBL_MAX;  // very large, but not infinity (so that infinite si values are discarded)
                    const int NO_INDEX = -99;  // meaningless cell index
                    int mq = NO_INDEX;

                    // plucker coords
                    const Plucker ray = Plucker(r(), k());

                    for (int face = 0; face < 4; face++)
                    {
                        std::array<Plucker, 3> edges;
                        Plucker::face(edges, _grid->_tetrahedra[_mr]->_vertices, face);

                        int leaveFace = -1;
                        for (int i = 0; i < 3; i++)
                        {
                            double prod = Plucker::dot(ray, edges[i]);
                            // how to ensure edges are aligned the same way?
                            // prod < 0 or prod > 0 we don't know
                        }
                    }

                    // loop over the list of neighbor indices
                    const vector<int>& mv = _grid->_cells[_mr]->neighbors();
                    int n = mv.size();
                    for (int i = 0; i < n; i++)
                    {
                        int mi = mv[i];

                        // declare the intersection distance for this neighbor (init to a value that will be rejected)
                        double si = 0;

                        // --- intersection with neighboring cell
                        if (mi >= 0)
                        {
                            // get the site position for this neighbor
                            Vec pi = _grid->_cells[mi]->position();

                            // calculate the (unnormalized) normal on the bisecting plane
                            Vec n = pi - pr;

                            // calculate the denominator of the intersection quotient
                            double ndotk = Vec::dot(n, k());

                            // if the denominator is negative the intersection distance is negative,
                            // so don't calculate it
                            if (ndotk > 0)
                            {
                                // calculate a point on the bisecting plane
                                Vec p = 0.5 * (pi + pr);

                                // calculate the intersection distance
                                si = Vec::dot(n, p - r()) / ndotk;
                            }
                        }

                        // --- intersection with domain wall
                        else
                        {
                            switch (mi)
                            {
                                case -1: si = (_grid->extent().xmin() - rx()) / kx(); break;
                                case -2: si = (_grid->extent().xmax() - rx()) / kx(); break;
                                case -3: si = (_grid->extent().ymin() - ry()) / ky(); break;
                                case -4: si = (_grid->extent().ymax() - ry()) / ky(); break;
                                case -5: si = (_grid->extent().zmin() - rz()) / kz(); break;
                                case -6: si = (_grid->extent().zmax() - rz()) / kz(); break;
                                default: throw FATALERROR("Invalid neighbor ID");
                            }
                        }

                        // remember the smallest nonnegative intersection point
                        if (si > 0 && si < sq)
                        {
                            sq = si;
                            mq = mi;
                        }
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
                        _mr = mq;

                        // if we're outside the domain, terminate the path after returning this path segment
                        if (_mr < 0) setState(State::Outside);
                        return true;
                    }
                }
            }

            case State::Outside:
            {
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
