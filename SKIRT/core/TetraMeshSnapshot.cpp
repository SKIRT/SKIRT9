/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "TetraMeshSnapshot.hpp"
#include "EntityCollection.hpp"
#include "FatalError.hpp"
#include "Log.hpp"
#include "MediumSystem.hpp"
#include "NR.hpp"
#include "Parallel.hpp"
#include "ParallelFactory.hpp"
#include "PathSegmentGenerator.hpp"
#include "ProcessManager.hpp"
#include "Random.hpp"
#include "SimulationItem.hpp"
#include "SiteListInterface.hpp"
#include "SpatialGridPath.hpp"
#include "SpatialGridPlotFile.hpp"
#include "StringUtils.hpp"
#include "Table.hpp"
#include "TextInFile.hpp"
#include "Units.hpp"
#include "tetgen.h"
#include <chrono>
#include <iostream>
#include <queue>
#include <unordered_set>
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
    Node* child(Vec bfr, const vector<Vec*>& points) const
    {
        return lessthan(bfr, *points[_m], _axis) ? _left : _right;
    }

    // returns the other child than the one that would be apropriate for the specified query point
    Node* otherChild(Vec bfr, const vector<Vec*>& points) const
    {
        return lessthan(bfr, *points[_m], _axis) ? _right : _left;
    }

    // returns the squared distance from the query point to the split plane
    double squaredDistanceToSplitPlane(Vec bfr, const vector<Vec*>& points) const
    {
        switch (_axis)
        {
            case 0:  // split on x
                return sqr(points[_m]->x() - bfr.x());
            case 1:  // split on y
                return sqr(points[_m]->y() - bfr.y());
            case 2:  // split on z
                return sqr(points[_m]->z() - bfr.z());
            default:  // this should never happen
                return 0;
        }
    }

    // returns the node in this subtree that represents the site nearest to the query point
    Node* nearest(Vec bfr, const vector<Vec*>& points)
    {
        // recursively descend the tree until a leaf node is reached, going left or right depending on
        // whether the specified point is less than or greater than the current node in the split dimension
        Node* current = this;
        while (Node* child = current->child(bfr, points)) current = child;

        // unwind the recursion, looking for the nearest node while climbing up
        Node* best = current;
        double bestSD = (*points[best->m()] - bfr).norm2();
        while (true)
        {
            // if the current node is closer than the current best, then it becomes the current best
            double currentSD = (*points[current->m()] - bfr).norm2();
            if (currentSD < bestSD)
            {
                best = current;
                bestSD = currentSD;
            }

            // if there could be points on the other side of the splitting plane for the current node
            // that are closer to the search point than the current best, then ...
            double splitSD = current->squaredDistanceToSplitPlane(bfr, points);
            if (splitSD < bestSD)
            {
                // move down the other branch of the tree from the current node looking for closer points,
                // following the same recursive process as the entire search
                Node* other = current->otherChild(bfr, points);
                if (other)
                {
                    Node* otherBest = other->nearest(bfr, points);
                    double otherBestSD = (*points[otherBest->m()] - bfr).norm2();
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
    for (auto cell : _vertices) delete cell;
    for (auto tetra : _tetrahedra) delete tetra;  // WIP FIX THIS
    for (auto tree : _blocktrees) delete tree;
}

////////////////////////////////////////////////////////////////////

void TetraMeshSnapshot::setExtent(const Box& extent)
{
    _extent = extent;
    _eps = 1e-12 * extent.widths().norm();
}

////////////////////////////////////////////////////////////////////

TetraMeshSnapshot::TetraMeshSnapshot(const TetraMeshSpatialGrid* grid, const Box& extent)
{
    setContext(grid);
    setExtent(extent);
    buildMesh(grid);
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

void TetraMeshSnapshot::buildMesh(const TetraMeshSpatialGrid* grid)
{
    tetgenio in, out;
    // tetgenio::facet* f;
    // tetgenio::polygon* p;

    in.firstnumber = 0;
    in.numberofpoints = 8;
    // in.numberofpoints += 4;
    in.pointlist = new REAL[in.numberofpoints * 3];

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
        // x, y the same but z = zmax
        in.pointlist[12 + i * 3 + 0] = in.pointlist[i * 3 + 0];
        in.pointlist[12 + i * 3 + 1] = in.pointlist[i * 3 + 1];
        in.pointlist[12 + i * 3 + 2] = _extent.zmax();
    }

    // for (int i = 0; i < 4; i++)
    // {
    //     in.pointlist[24 + i * 3 + 0] = in.pointlist[i * 3 + 0] + _extent.xwidth() * 0.2;
    //     in.pointlist[24 + i * 3 + 1] = in.pointlist[i * 3 + 1] + _extent.ywidth() * 0.2;
    //     in.pointlist[24 + i * 3 + 2] = in.pointlist[i * 3 + 2] + _extent.zwidth() * 0.2;
    // }

    /* psc */
    // in.numberofpointattributes = 1;
    // in.pointattributelist = new REAL[in.numberofpoints];
    // for (int i = 0; i < 8; i++)
    // {
    //     in.pointattributelist[i] = 1.;
    // }
    // for (int i = 8; i < in.numberofpoints; i++)
    // {
    //     Vec pos = random()->position(_extent);
    //     in.pointlist[i * 3 + 0] = pos.x();
    //     in.pointlist[i * 3 + 1] = pos.y();
    //     in.pointlist[i * 3 + 2] = pos.z();
    // }
    // for (int i = 8; i < in.numberofpoints; i++)
    // {
    //     in.pointattributelist[i] = 1.;
    // }
    // in.pointattributelist[8] = 50.;

    in.numberoffacets = 6;
    // in.numberoffacets += 1;
    in.facetlist = new tetgenio::facet[in.numberoffacets];
    addFacet(&in.facetlist[0], {0, 1, 2, 3});  // Facet 1. bottom
    addFacet(&in.facetlist[1], {4, 5, 6, 7});  // Facet 2. top
    addFacet(&in.facetlist[2], {0, 4, 5, 1});  // Facet 3. front
    addFacet(&in.facetlist[3], {1, 5, 6, 2});  // Facet 4. right
    addFacet(&in.facetlist[4], {2, 6, 7, 3});  // Facet 5. back
    addFacet(&in.facetlist[5], {3, 7, 4, 0});  // Facet 6. left
    // addFacet(&in.facetlist[6], {8, 9, 10, 11});  // Facet 6. left

    tetgenbehavior behavior;
    behavior.plc = 1;          // -p PLC
    behavior.quality = 1;      // -q quality mesh
    behavior.fixedvolume = 1;  // -a max volume
    behavior.neighout = 2;     // -nn neighbors and edges?
    behavior.zeroindex = 1;    // -z zero index
    // behavior.edgesout = 1;     // -e edges
    // behavior.weighted = 1;     // -w weighted

    // parameters
    // behavior.minratio = 5.0;   // -q quality
    behavior.maxvolume = 0.05 * _extent.volume();  // -a max volume
    // behavior.weighted_param = ;
    // behavior.mindihedral = 5.0;  // -q/ minimal angle

    in.tetunsuitable = [grid](double* pa, double* pb, double* pc, double* pd, double vol) {
        return grid->tetUnsuitable(pa, pb, pc, pd, vol);
    };

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

    _tetrahedra.resize(numTetra);
    for (int i = 0; i < numTetra; i++)
    {
        std::array<Vec*, 4> vertices;
        std::array<Face, 4> neighbors;
        for (int c = 0; c < 4; c++)
        {
            vertices[c] = _vertices[out.tetrahedronlist[4 * i + c]];
            int ntetra = out.neighborlist[4 * i + c];

            // find which face is shared with neighbor
            int nface;
            for (int cn = 0; cn < 4; cn++)
            {
                if (out.neighborlist[4 * ntetra + cn] == i)
                {
                    nface = cn;
                    break;
                }
            }
            neighbors[c] = Face(ntetra, nface);
        }

        _tetrahedra[i] = new Tetra(vertices, neighbors);

        _centroids.push_back(&_tetrahedra[i]->_centroid);
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
    _rhov.resize(numVertices);
    Array Mv(numVertices);

    // get the maximum temperature, or zero of there is none
    double maxT = useTemperatureCutoff() ? maxTemperature() : 0.;

    // initialize statistics
    double totalOriginalMass = 0;
    double totalMetallicMass = 0;
    double totalEffectiveMass = 0;

    // loop over all sites/cells
    int numIgnored = 0;
    for (int m = 0; m != numVertices; ++m)
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
    if (numVertices) NR::cdf(_cumrhov, Mv);
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
            return m1 != m2 && lessthan(*_centroids[m1], *_centroids[m2], depth % 3);
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
    if (!numTetra) return;

    log()->info("Building data structures to accelerate searching the tetrahedralization");

    // -------------  block lists  -------------
    _nb = max(3, min(250, static_cast<int>(cbrt(numTetra))));
    // _nb = 1;
    _nb2 = _nb * _nb;
    _nb3 = _nb * _nb * _nb;

    // initialize a vector of nb * nb * nb lists
    _blocklists.resize(_nb3);

    // we add the tetrahedra to all blocks they potentially overlap with
    // this will slow down the search tree but if no search tree is present
    // we can simply loop over all tetrahedra inside the block
    int i1, j1, k1, i2, j2, k2;
    for (int c = 0; c != numTetra; ++c)
    {
        // _extent.cellIndices(i1, j1, k1, *_centroids[c], _nb, _nb, _nb);
        // _blocklists[i1 * _nb2 + j1 * _nb + k1].push_back(c);

        _extent.cellIndices(i1, j1, k1, _tetrahedra[c]->rmin() - Vec(_eps, _eps, _eps), _nb, _nb, _nb);
        _extent.cellIndices(i2, j2, k2, _tetrahedra[c]->rmax() + Vec(_eps, _eps, _eps), _nb, _nb, _nb);
        for (int i = i1; i <= i2; i++)
            for (int j = j1; j <= j2; j++)
                for (int k = k1; k <= k2; k++) _blocklists[i * _nb2 + j * _nb + k].push_back(c);
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
    log()->info("  Number of sites: " + std::to_string(numVertices));

    // abort if there are no cells
    if (!numVertices) return;

    // construct a single search tree on the site locations of all cells
    log()->info("Building data structure to accelerate searching " + std::to_string(numVertices) + " Voronoi sites");
    _blocktrees.resize(1);
    vector<int> ids(numVertices);
    for (int m = 0; m != numVertices; ++m) ids[m] = m;
    _blocktrees[0] = buildTree(ids.begin(), ids.end(), 0);
}

////////////////////////////////////////////////////////////////////

void TetraMeshSnapshot::writeGridPlotFiles(const SimulationItem* probe) const
{
    ///////////////// TEMP
    std::ofstream outputFile("data/tetrahedra.txt");
    for (size_t i = 0; i < _tetrahedra.size(); i++)
    {
        const Tetra* tetra = _tetrahedra[i];
        for (size_t l = 0; l < 4; l++)
        {
            const Vec* r = tetra->_vertices[l];
            outputFile << r->x() << ", " << r->y() << ", " << r->z() << "\n";
        }
    }
    outputFile.close();
    /////////////////

    // create the plot files
    SpatialGridPlotFile plotxy(probe, probe->itemName() + "_grid_xy");
    SpatialGridPlotFile plotxz(probe, probe->itemName() + "_grid_xz");
    SpatialGridPlotFile plotyz(probe, probe->itemName() + "_grid_yz");
    SpatialGridPlotFile plotxyz(probe, probe->itemName() + "_grid_xyz");

    // for each site, compute the corresponding cell and output its edges
    log()->info("Writing plot files for tetrahedralization with " + std::to_string(numTetra) + " tetrahedra");
    log()->infoSetElapsed(numTetra);
    int numDone = 0;
    for (int i = 0; i < numTetra; i++)
    {
        const Tetra* tetra = _tetrahedra[i];
        vector<double> coords;
        coords.reserve(12);
        vector<int> indices;
        indices.reserve(16);

        for (int v = 0; v < 4; v++)
        {
            const Vec* vertex = tetra->_vertices[v];
            coords.push_back(vertex->x());
            coords.push_back(vertex->y());
            coords.push_back(vertex->z());

            // get vertices of opposite face
            std::array<int, 3> faceIndices = tetra->clockwiseVertices(v);
            indices.push_back(3);  // amount of vertices per face
            indices.push_back(faceIndices[0]);
            indices.push_back(faceIndices[1]);
            indices.push_back(faceIndices[2]);
        }

        if (tetra->zmin() <= 0 && tetra->zmax() >= 0) plotxy.writePolyhedron(coords, indices);
        if (tetra->ymin() <= 0 && tetra->ymax() >= 0) plotxz.writePolyhedron(coords, indices);
        if (tetra->xmin() <= 0 && tetra->xmax() >= 0) plotyz.writePolyhedron(coords, indices);
        if (i <= 1000)
            plotxyz.writePolyhedron(coords, indices);  // like VoronoiMeshSnapshot, but why even write at all?

        // log message if the minimum time has elapsed
        numDone++;
        if (numDone % 2000 == 0) log()->infoIfElapsed("Computed tetrehedra: ", 2000);
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
    return Position(_tetrahedra[m]->_centroid);
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
    double s = random()->uniform();
    double t = random()->uniform();
    double u = random()->uniform();
    return _tetrahedra[m]->generatePosition(s, t, u);
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
    if (!_extent.contains(bfr)) return -1;

    // determine the block in which the point falls
    int i, j, k;
    _extent.cellIndices(i, j, k, bfr, _nb, _nb, _nb);
    int b = i * _nb2 + j * _nb + k;

    // look for the closest centroid in this block using the search tree
    Node* tree = _blocktrees[b];
    if (tree)
    {
        // use a breadth-first search over the neighboring tetrahedra
        int root = tree->nearest(bfr, _centroids)->m();
        std::queue<int> queue;
        std::unordered_set<int> explored;
        queue.push(root);

        int t;
        while (queue.size() > 0)
        {
            t = queue.front();
            queue.pop();
            const Tetra* tetra = _tetrahedra[t];

            if (tetra->inside(bfr)) return t;
            explored.insert(t);

            for (const Face& face : tetra->_faces)
            {
                int n = face._ntetra;
                if (n != -1 && explored.find(n) == explored.end())  // if not already explored
                    queue.push(n);
            }
        }
        log()->error("cellIndex failed to find the tetrahedron");  // change to warning?
    }
    else
    {
        // if there is no search tree, simply loop over all tetrahedra in the block
        for (int t : _blocklists[b])
        {
            if (_tetrahedra[t]->inside(bfr)) return t;
        }
    }

    return -1;
}

////////////////////////////////////////////////////////////////////

const Array& TetraMeshSnapshot::properties(int /*m*/) const
{
    // return _sites[m]->properties();
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
    int _mr = -1;
    int _enteringFace = -1;
    int _leavingFace = -1;

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

    void write(double ds, int face)
    {
        out << "photon=" << _mr << "," << face << "\n";
        out << "ds=" << ds << "\n";
        out << "r=" << r().x() << "," << r().y() << "," << r().z() << "\n";
        out << "k=" << k().x() << "," << k().y() << "," << k().z() << std::endl;
    }
#endif
    bool next() override
    {
        if (state() == State::Unknown)
        {
            // try moving the photon packet inside the grid; if this is impossible, return an empty path
            // this also changes the state()
            if (!moveInside(_grid->extent(), _grid->_eps)) return false;

            // get the index of the cell containing the current position
            _mr = _grid->cellIndex(r());

            // if the photon packet started outside the grid, return the corresponding nonzero-length segment;
            // otherwise fall through to determine the first actual segment
            // if (ds() > 0.) return true;  // is this really needed?

            // find the entering face, this will do a few full plucker products for the first traversal step
            if (state() == State::Inside)
            {

                const Tetra* tetra = _grid->_tetrahedra[_mr];
                // Plücker coordinates of the ray
                const Vec U = k();
                const Vec V = Vec::cross(U, r());

                for (int face = 0; face < 4; face++)
                {
                    bool entering = true;
                    std::array<int, 3> cv = tetra->clockwiseVertices(face);

                    for (int i = 0; i < 3; i++)
                    {
                        // edges: 12, 20, 01
                        // verts:  0,  1,  2
                        int t1 = cv[(i + 1) % 3];
                        int t2 = cv[(i + 2) % 3];
                        Vec edge_U = tetra->getEdge(t1, t2);
                        Vec edge_V = Vec::cross(edge_U, *tetra->_vertices[t1]);

                        // naive implementation for Plücker product with all 3 edges per face
                        double prod = Vec::dot(U, edge_V) + Vec::dot(V, edge_U);

                        // all products must be non-negative for it to be an entering face
                        if (prod < 0)
                        {
                            entering = false;
                            break;
                        }
                    }

                    if (entering)
                    {
                        _enteringFace = face;
                        break;
                    }
                }
            }
#ifdef WRITE
            startwriting();
            write(0, _mr);
#endif
        }

        // intentionally falls through
        if (state() == State::Inside)
        {
            // loop in case no exit point was found (which should happen only rarely)
            while (true)
            {
                const Tetra* tetra = _grid->_tetrahedra[_mr];

                // the translated Plücker moment in the local coordinate system
                const Vec moment = Vec::cross(k(), r() - *tetra->_vertices[_enteringFace]);
                std::array<int, 3> cv = Tetra::clockwiseVertices(_enteringFace);

                // 2 step decision tree
                int clockwise0 = Vec::dot(moment, tetra->getEdge(cv[0], _enteringFace)) <= 0; // problem here is prod = 0
                // if clockwise move clockwise
                // (0+1)%3=1 is c and (0-1)%3=2 is cc
                int i = clockwise0 ? 1 : 2;
                int clockwisei = Vec::dot(moment, tetra->getEdge(cv[i], _enteringFace)) < 0;
                // c0 and ci form binary % 3
                //  0 0 -> 0
                //  0 1 -> 1
                //  1 0 -> 2
                //  1 1 -> 0

                _leavingFace = cv[((clockwise0 << 1) | clockwisei) % 3];

                // calculate ds from exit face
                const Vec& v0 = *tetra->_vertices[0];
                Vec e01 = tetra->getEdge(0, 1);
                Vec e02 = tetra->getEdge(0, 2);
                Vec n = Vec::cross(e01, e02);

                double ds = Vec::dot(n, v0 - r()) / Vec::dot(n, k());

                // we could assert if r+exit and r+sq*k are the same

                propagater(ds + _grid->_eps);
                setSegment(_mr, ds);
#ifdef WRITE
                write(ds, _leavingFace);
#endif

                _mr = tetra->_faces[_leavingFace]._ntetra;

                // set enteringFace only if there is a neighboring cell
                if (_mr < 0)
                {
                    // _enteringFace = -1;
                    setState(State::Outside);
                }
                else
                {
                    _enteringFace = tetra->_faces[_leavingFace]._nface;
                }

                return true;
            }
        }

        if (state() == State::Outside)
        {}
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