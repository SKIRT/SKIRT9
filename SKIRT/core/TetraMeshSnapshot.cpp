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
#include "tetgen.h"
#include <iostream>
#include <set>
#include "container.hh"

////////////////////////////////////////////////////////////////////

namespace
{
    // classes used for serializing/deserializing Voronoi cell geometry when communicating the results of
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
class TetraMeshSnapshot::Site
{
public:
    Vec _r;                  // site position
    vector<int> _neighbors;  // list of neighbor indices in _cells vector
    Array _properties;       // user-defined properties, if any

public:
    // constructor stores the specified site position; the other data members are set to zero or empty
    Site(Vec r) : _r(r) {}

    // constructor derives the site position from the first three property values and stores the user properties;
    // the other data members are set to zero or empty
    Site(const Array& prop) : _r(prop[0], prop[1], prop[2]), _properties{prop} {}  // WIP

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

    // writes the Voronoi cell geometry to the serialized data buffer, preceded by the specified cell index,
    // if the cell geometry has been calculated for this cell; otherwise does nothing
    void writeGeometryIfPresent(SerializedWrite& wdata, int m)
    {
        if (!_neighbors.empty())
        {
            wdata.write(m);
            wdata.write(_neighbors);
        }
    }

    // reads the Voronoi cell geometry from the serialized data buffer
    void readGeometry(SerializedRead& rdata) { rdata.read(_neighbors); }
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

TetraMeshSnapshot::Edge::Edge(int i1, int i2, const Vec* v1, const Vec* v2) : Plucker(*v1, *v2 - *v1), i1(i1), i2(i2) {}

////////////////////////////////////////////////////////////////////

TetraMeshSnapshot::Tetra::Tetra(const std::array<Vec*, 4>& vertices, const std::array<int, 4>& indices,
                                const std::array<int, 4>& neighbors, const std::array<Edge*, 6>& edges)
    : _vertices(vertices), _indices(indices), _neighbors(neighbors), _edges(edges)
{
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

    const Vec e01 = *_vertices[1] - *_vertices[0];
    const Vec e02 = *_vertices[2] - *_vertices[0];
    const Vec e03 = *_vertices[3] - *_vertices[0];
    // this convention makes edges go clockwise around leaving rays from inside the tetrahedron
    // so their plucker products are all positive if the ray leaves
    double orientation = Vec::dot(Vec::cross(e01, e02), e03);
    if (orientation < 0)
    {
        std::cout << "ORIENTATION SWITCHED!!!!!!!!!" << std::endl;
        // swap last 2, this means first 2 indices can be ordered i < j
        std::swap(_vertices[2], _vertices[3]);
        std::swap(_neighbors[2], _neighbors[3]);
        std::swap(_indices[2], _indices[3]);
    }
}

////////////////////////////////////////////////////////////////////

double TetraMeshSnapshot::Tetra::getProd(const Plucker& ray, int t1, int t2) const
{
    int e = (std::min(t1, t2) == 0) ? std::max(t1, t2) - 1 : t1 + t2;
    Edge* edge = _edges[e];

    // not the same order -> *-1
    // beginning of t1 == beginning of edge (= same order)
    return (_indices[t1] == edge->i1 ? 1 : -1) * Plucker::dot(ray, *edge);
}

////////////////////////////////////////////////////////////////////

bool TetraMeshSnapshot::Tetra::intersects(std::array<double, 3>& barycoords, const Plucker& ray, int face,
                                          bool leaving) const
{
    std::array<int, 3> t = counterclockwiseVertices(face);

    double sum = 0;
    for (int i = 0; i < 3; i++)
    {
        // edges: 12, 20, 01
        // verts:  0,  1,  2
        double prod = getProd(ray, t[(i + 1) % 3], t[(i + 2) % 3]);
        if (leaving != (prod >= 0)) return false;  // change this so for both leavig and entering prod=0 works
        barycoords[i] = prod;
        sum += prod;
    }
    for (int i = 0; i < 3; i++) barycoords[i] /= sum;

    return true;
}

////////////////////////////////////////////////////////////////////

bool TetraMeshSnapshot::Tetra::inside(const Position& bfr) const
{
    if (!Box::contains(bfr)) return false;

    /*
    face: normals for which the other vertex has a positive dot product with
    3:*02 x 01*| 10 x 12 | 21 x 20
    2: 13 x 10 |*01 x 03*| 30 x 31
    1: 20 x 23 | 32 x 30 |*03 x 02*
    0: 31 x 32 | 23 x 21 |*12 x 13* // last one doesn't matter
    */

    // optimized version (probably not that much better)
    Vec e0p = bfr - *_vertices[0];
    Vec e02 = *_vertices[2] - *_vertices[0];
    Vec e01 = *_vertices[1] - *_vertices[0];
    if (Vec::dot(Vec::cross(e02, e01), e0p) > 0)  // 02 x 01
        return false;

    Vec e03 = *_vertices[3] - *_vertices[0];
    if (Vec::dot(Vec::cross(e01, e03), e0p) > 0)  // 01 x 03
        return false;

    if (Vec::dot(Vec::cross(e03, e02), e0p) > 0)  // 03 x 02
        return false;

    Vec e1p = bfr - *_vertices[1];
    Vec e12 = *_vertices[2] - *_vertices[1];
    Vec e13 = *_vertices[3] - *_vertices[1];
    return Vec::dot(Vec::cross(e12, e13), e1p) < 0;  // 12 x 13

    // checks 3 edges too many but very simple
    // for (int face = 0; face < 4; face++)
    // {
    //     std::array<int, 3> t = clockwiseVertices(face);
    //     Vec& v0 = *_vertices[t[0]];
    //     Vec& clock = *_vertices[t[1]];
    //     Vec& counter = *_vertices[t[2]];
    //     Vec normal = Vec::cross(counter - v0, clock - v0);
    //     if (Vec::dot(normal, bfr - v0) < 0)  // is pos on the same side as v3
    //     {
    //         return false;
    //     }
    // }
    // return true;
}

////////////////////////////////////////////////////////////////////

Vec TetraMeshSnapshot::Tetra::calcExit(const std::array<double, 3>& barycoords, int face) const
{
    std::array<int, 3> t = Tetra::counterclockwiseVertices(face);
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

std::array<int, 3> TetraMeshSnapshot::Tetra::counterclockwiseVertices(int face)
{
    std::array<int, 3> t = {(face + 1) % 4, (face + 2) % 4, (face + 3) % 4};
    // if face is even we should swap two edges
    if (face % 2 == 0) std::swap(t[0], t[2]);
    return t;
}

////////////////////////////////////////////////////////////////////

namespace
{
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

    void addFacet(tetgenio::facet* f, std::array<int, 4> vertices)
    {
        f->numberofpolygons = 1;
        f->polygonlist = new tetgenio::polygon[f->numberofpolygons];
        f->numberofholes = 0;
        f->holelist = NULL;
        tetgenio::polygon* p = &f->polygonlist[0];
        p->numberofvertices = 4;
        p->vertexlist = new int[p->numberofvertices];
        p->vertexlist[0] = vertices[0];
        p->vertexlist[1] = vertices[1];
        p->vertexlist[2] = vertices[2];
        p->vertexlist[3] = vertices[3];
    }
}

////////////////////////////////////////////////////////////////////

class TetraMeshSnapshot::Node
{
private:
    int _m;        // index in _cells to the site defining the split at this node
    int _axis;     // split axis for this node (0,1,2)
    Node* _up;     // ptr to the parent node
    Node* _left;   // ptr to the left child node
    Node* _right;  // ptr to the right child node
    std::vector<int> tetra;

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
    Node* child(Vec bfr, const vector<Site*>& sites) const
    {
        return lessthan(bfr, sites[_m]->position(), _axis) ? _left : _right;
    }

    // returns the other child than the one that would be apropriate for the specified query point
    Node* otherChild(Vec bfr, const vector<Site*>& sites) const
    {
        return lessthan(bfr, sites[_m]->position(), _axis) ? _right : _left;
    }

    // returns the squared distance from the query point to the split plane
    double squaredDistanceToSplitPlane(Vec bfr, const vector<Site*>& sites) const
    {
        switch (_axis)
        {
            case 0:  // split on x
                return sqr(sites[_m]->position().x() - bfr.x());
            case 1:  // split on y
                return sqr(sites[_m]->position().y() - bfr.y());
            case 2:  // split on z
                return sqr(sites[_m]->position().z() - bfr.z());
            default:  // this should never happen
                return 0;
        }
    }

    // returns the node in this subtree that represents the site nearest to the query point
    Node* nearest(Vec bfr, const vector<Site*>& sites)
    {
        // recursively descend the tree until a leaf node is reached, going left or right depending on
        // whether the specified point is less than or greater than the current node in the split dimension
        Node* current = this;
        while (Node* child = current->child(bfr, sites)) current = child;

        // unwind the recursion, looking for the nearest node while climbing up
        Node* best = current;
        double bestSD = sites[best->m()]->squaredDistanceTo(bfr);
        while (true)
        {
            // if the current node is closer than the current best, then it becomes the current best
            double currentSD = sites[current->m()]->squaredDistanceTo(bfr);
            if (currentSD < bestSD)
            {
                best = current;
                bestSD = currentSD;
            }

            // if there could be points on the other side of the splitting plane for the current node
            // that are closer to the search point than the current best, then ...
            double splitSD = current->squaredDistanceToSplitPlane(bfr, sites);
            if (splitSD < bestSD)
            {
                // move down the other branch of the tree from the current node looking for closer points,
                // following the same recursive process as the entire search
                Node* other = current->otherChild(bfr, sites);
                if (other)
                {
                    Node* otherBest = other->nearest(bfr, sites);
                    double otherBestSD = sites[otherBest->m()]->squaredDistanceTo(bfr);
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
    for (auto cell : _sites) delete cell;
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
    while (in.readRow(coords)) _sites.push_back(new Site(Vec(coords[0], coords[1], coords[2])));
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
    _sites.resize(n);
    for (int m = 0; m != n; ++m) _sites[m] = new Site(sli->sitePosition(m));

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
    _sites.resize(n);
    for (int m = 0; m != n; ++m) _sites[m] = new Site(sites[m]);

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
    tetgenio in, out;
    tetgenio::facet* f;
    tetgenio::polygon* p;

    in.firstnumber = 0;
    in.numberofpoints = 8;
    in.pointlist = new REAL[in.numberofpoints * 3];

    // _extent = Box(-1, -1, -1, 1, 1, 1);

    // bottom half (zmin)
    in.pointlist[0] = _extent.xmin();
    in.pointlist[1] = _extent.ymin();
    in.pointlist[2] = _extent.zmin();

    in.pointlist[3] = _extent.xmax();
    in.pointlist[4] = _extent.xmin();
    in.pointlist[5] = _extent.zmin();

    in.pointlist[6] = _extent.xmax();
    in.pointlist[7] = _extent.xmax();
    in.pointlist[8] = _extent.zmin();

    in.pointlist[9] = _extent.xmin();
    in.pointlist[10] = _extent.ymax();
    in.pointlist[11] = _extent.zmin();

    // top half (zmax)
    for (int i = 0; i < 4; i++)
    {
        in.pointlist[12 + i * 3 + 0] = in.pointlist[i * 3 + 0];
        in.pointlist[12 + i * 3 + 1] = in.pointlist[i * 3 + 1];
        in.pointlist[12 + i * 3 + 2] = _extent.zmax();
    }

    in.numberoffacets = 6;
    in.facetlist = new tetgenio::facet[in.numberoffacets];
    addFacet(&in.facetlist[0], {0, 1, 2, 3});  // Facet 1. bottom
    addFacet(&in.facetlist[1], {4, 5, 6, 7});  // Facet 2. top
    addFacet(&in.facetlist[2], {0, 4, 5, 1});  // Facet 3. front
    addFacet(&in.facetlist[3], {1, 5, 6, 2});  // Facet 4. right
    addFacet(&in.facetlist[4], {2, 6, 7, 3});  // Facet 5. back
    addFacet(&in.facetlist[5], {3, 7, 4, 0});  // Facet 6. left

    tetgenbehavior behavior;
    behavior.plc = 1;        // -p PLC
    behavior.quality = 1;    // -q quality mesh
    behavior.neighout = 2;   // -nn neighbors and edges?
    behavior.zeroindex = 1;  // -z zero index
    behavior.edgesout = 1;   // -e edges

    // parameters
    behavior.minratio = 2.0;   // -q quality
    behavior.maxvolume = 0.9;  // -a max volume
    // behavior.mindihedral = 5.0;  // -q/ minimal angle

    tetrahedralize(&behavior, &in, &out);
    numTetra = out.numberoftetrahedra;
    numEdges = out.numberofedges;
    numVertices = out.numberofpoints;

    _vertices.resize(numVertices);
    for (int i = 0; i < numVertices; i++)
    {
        double x = out.pointlist[3 * i + 0];
        double y = out.pointlist[3 * i + 1];
        double z = out.pointlist[3 * i + 2];

        _vertices[i] = new Vec(x, y, z);
    }

    _edges.resize(numEdges);
    for (int i = 0; i < numEdges; i++)
    {
        int v1 = out.edgelist[2 * i];
        int v2 = out.edgelist[2 * i + 1];
        _edges[i] = new Edge(v1, v2, _vertices[v1], _vertices[v2]);
    }

    _tetrahedra.resize(numTetra);
    for (int i = 0; i < numTetra; i++)
    {
        std::array<Vec*, 4> vertices;
        std::array<int, 4> indices;
        std::array<int, 4> neighbors;
        std::array<Edge*, 6> edges;
        for (int c = 0; c < 4; c++)
        {
            indices[c] = out.tetrahedronlist[4 * i + c];
            vertices[c] = _vertices[indices[c]];
            neighbors[c] = out.neighborlist[4 * i + c];
        }
        for (int e = 0; e < 6; e++)
        {
            int ei = out.tet2edgelist[6 * i + e];

            Edge* edge = _edges[ei];
            auto t1 = std::find(indices.begin(), indices.end(), edge->i1) - indices.begin();
            auto t2 = std::find(indices.begin(), indices.end(), edge->i2) - indices.begin();

            // tetgen edge order: 23 03 01 12 13 02
            // static constexpr int tetgen_order[12] = {2, 3, 0, 3, 0, 1, 1, 2, 1, 3, 0, 2};
            // int t1 = tetgen_order[2 * e];      // 2 0 0 1 1 0
            // int t2 = tetgen_order[2 * e + 1];  // 3 3 1 2 3 2

            if (t1 > t2) std::swap(t1, t2);
            edges[(t1 == 0) ? t2 - 1 : t1 + t2] = _edges[ei];
        }
        _tetrahedra[i] = new Tetra(vertices, indices, neighbors, edges);
    }

    log()->info("number of vertices " + std::to_string(numVertices));
    log()->info("number of edges " + std::to_string(numEdges));
    log()->info("number of tetrahedra " + std::to_string(numTetra));
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

TetraMeshSnapshot::Node* TetraMeshSnapshot::buildTree(vector<int>::iterator first, vector<int>::iterator last,
                                                      int depth) const
{
    auto length = last - first;
    if (length > 0)
    {
        auto median = length >> 1;
        std::nth_element(first, first + median, last, [this, depth](int m1, int m2) {
            return m1 != m2 && lessthan(_sites[m1]->position(), _sites[m2]->position(), depth % 3);
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
    int numTetra = _sites.size();
    if (!numTetra) return;

    log()->info("Building data structures to accelerate searching the Tetra tesselation");

    // -------------  block lists  -------------

    // initialize a vector of nb x nb x nb lists, each containing the cells overlapping a certain block in the domain
    _blocklists.resize(_nb3);

    // add the cell object to the lists for all blocks it may overlap
    int i1, j1, k1, i2, j2, k2;
    for (int m = 0; m != numTetra; ++m)
    {
        auto& vert = _tetrahedra[m]->_vertices;
        auto vert_minmax =
            std::minmax_element(vert.begin(), vert.end(), [](Vec* a, Vec* b) { return a->norm2() < b->norm2(); });

        _extent.cellIndices(i1, j1, k1, **vert_minmax.first - Vec(_eps, _eps, _eps), _nb, _nb, _nb);
        _extent.cellIndices(i2, j2, k2, **vert_minmax.second + Vec(_eps, _eps, _eps), _nb, _nb, _nb);
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
    int numCells = _sites.size();
    log()->info("  Number of sites: " + std::to_string(numCells));

    // abort if there are no cells
    if (!numCells) return;

    // construct a single search tree on the site locations of all cells
    log()->info("Building data structure to accelerate searching " + std::to_string(numCells) + " Voronoi sites");
    _blocktrees.resize(1);
    vector<int> ids(numCells);
    for (int m = 0; m != numCells; ++m) ids[m] = m;
    _blocktrees[0] = buildTree(ids.begin(), ids.end(), 0);
}

////////////////////////////////////////////////////////////////////

void TetraMeshSnapshot::writeGridPlotFiles(const SimulationItem* /*probe*/) const
{
    std::ofstream outputFile("data/tetrahedra.txt");
    // outputFile << "vertices=" << _vertices.size() << "\n";
    // for (size_t i = 0; i < _vertices.size(); i++)
    // {
    //     outputFile << "vertex=" << i << "\n";
    //     Vec& r = *_vertices[i];
    //     outputFile << r.x() << ", " << r.y() << ", " << r.z() << "\n";
    //     outputFile << "\n";
    // }

    for (size_t i = 0; i < _tetrahedra.size(); i++)
    {
        const Tetra* tetra = _tetrahedra[i];
        for (size_t l = 0; l < 4; l++)
        {
            const Vec* r = tetra->_vertices[l];
            outputFile << r->x() << ", " << r->y() << ", " << r->z() << "\n";
        }
        // outputFile << "neighbors=";
        // for (size_t j = 0; j < 4; j++)
        // {
        //     outputFile << " " << tetra->_neighbors[j];
        //     if (j != 3) outputFile << ",";
        // }
        // outputFile << "\n";
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
    // make sure the position is inside the domain
    // if (!_extent.contains(bfr)) return -1;

    // determine the block in which the point falls
    // if we didn't build a Voronoi mesh, the search tree is always in the first "block"
    // int i, j, k;
    // _extent.cellIndices(i, j, k, bfr, _nb, _nb, _nb);
    // int b = i * _nb2 + j * _nb + k;

    // look for the closest site in this block, using the search tree if there is one
    // Node* tree = _blocktrees[b];
    // int m = tree->nearest(bfr, _sites)->m();
    // find all edges that connect to m
    // Site* site = _sites[m];

    // I can retreive all edges by looking at neighbors of Site
    // i want tetra though so I'll build hashmap for Tetra instead of vector

    // if there is no search tree, simply loop over the index list
    // maybe use a k-d tree here to find nearest tetrahedra
    for (int i = 0; i < numTetra; i++)
    {
        const Tetra* tetra = _tetrahedra[i];
        if (tetra->Tetra::inside(bfr)) return i;
    }
    return -1;
}

////////////////////////////////////////////////////////////////////

const Array& TetraMeshSnapshot::properties(int m) const
{
    return _sites[m]->properties();
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

#define WRITE

#ifdef WRITE
    std::ofstream out;

    void startwriting()
    {
        if (!out.is_open()) out.open("data/photon.txt");
    }

    void stopwriting()
    {
        out.close();
    }

    void write(const Vec& exit, double ds, int face)
    {
        out << "photon=" << _mr << "," << face << "\n";
        out << "ds=" << ds << "\n";
        out << "r=" << r().x() << "," << r().y() << "," << r().z() << "\n";
        out << "exit=" << exit.x() << "," << exit.y() << "," << exit.z() << "\n";
        out << "k=" << k().x() << "," << k().y() << "," << k().z() << std::endl;
    }
#endif
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

#ifdef WRITE
            startwriting();
            write(r(), 0, _mr);
#endif
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
                if (enteringFace == -1)  // start ray traversal inside
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
                    std::array<int, 3> t = Tetra::counterclockwiseVertices(enteringFace);

                    // 2 step decision tree
                    prods[0] = tetra->getProd(ray, t[0], enteringFace);
                    bool clockwise0 = prods[0] > 0;
                    int i = clockwise0 ? 1 : 2;  // if (counter)clockwise move (counter)clockwise
                    prods[i] = tetra->getProd(ray, t[i], enteringFace);

                    // if 2 clockwise: face=t0
                    // if 2 counter: face=t0
                    // if clockwise then counter: face=t2
                    // if counter then clockwise: face=t1

                    if (clockwise0 == (prods[i] > 0))
                        leavingFace = t[0];
                    else if (clockwise0)
                        leavingFace = t[2];
                    else
                        leavingFace = t[1];

                    // get prods (this calculates prod[i] again but is much cleaner)
                    tetra->intersects(prods, ray, leavingFace, true);
                }

                // if no exit point was found, advance the current point by a small distance,
                // recalculate the cell index, and return to the start of the loop
                if (leavingFace == -1)
                {
                    propagater(_grid->_eps);
                    _mr = _grid->cellIndex(r());
                    enteringFace = -1;

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
#ifdef WRITE
                    write(exit, ds, leavingFace);
#endif
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
                        setState(State::Outside);
                    }
                    // set new cell
                    _mr = next_mr;

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
#ifdef WRITE
                write(r(), ds, -1);
                stopwriting();
#endif
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
#ifdef WRITE
                                write(exit, ds, enteringFace);
#endif

                                wasInside = true;
                                setState(State::Inside);
                                return true;
                            }
                        }
                    }
                }
            }
        }
#ifdef WRITE
        stopwriting();
#endif
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
