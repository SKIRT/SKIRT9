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

    template<typename T> int sgn(T a)
    {
        return (a > T(0)) - (a < T(0));
    }

    int findEnteringFace(const Tetra* tetra, const Vec& pos, const Direction& dir)
    {
        int enteringFace = -1;
        // clockwise and cclockwise adjacent faces when checking edge v1->v2
        constexpr int etable[6][2] = {{3, 2}, {1, 3}, {2, 1}, {3, 0}, {0, 2}, {1, 0}};
        // try all 6 edges because of rare edge cases where ray is inside edge
        // having only 1 non-zero Plücker product
        int e = 0;
        for (int v1 = 0; v1 < 3; v1++)
        {
            for (int v2 = v1 + 1; v2 < 4; v2++)
            {

                Vec moment12 = Vec::cross(dir, pos - *tetra->_vertices[v1]);
                double prod12 = Vec::dot(moment12, tetra->getEdge(v1, v2));
                if (prod12 != 0.)
                {
                    enteringFace = prod12 < 0 ? etable[e][0] : etable[e][1];
                    break;
                }
                e++;
            }
            if (enteringFace != -1) break;
        }
        return enteringFace;
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

void TetraMeshSnapshot::readAndClose()
{
    // // read the site info into memory
    // Array prop;
    // while (infile()->readRow(prop))
    // {
    //     _vertices.push_back(new Vec(prop[0], prop[1], prop[2]));
    // }

    // // close the file
    // Snapshot::readAndClose();

    // // if we are allowed to build a Voronoi mesh
    // // calculate the Voronoi cells
    // buildMesh(?, false);

    // // if a mass density policy has been set, calculate masses and densities and build the search data structure
    // if (hasMassDensityPolicy()) calculateDensityAndMass();
    // if (hasMassDensityPolicy() || needGetEntities()) buildSearchPerBlock();
}

////////////////////////////////////////////////////////////////////

void TetraMeshSnapshot::setExtent(const Box& extent)
{
    _extent = extent;
    _eps = 1e-12 * extent.widths().norm();
}

////////////////////////////////////////////////////////////////////

// TetraMeshSnapshot::TetraMeshSnapshot(const SimulationItem* item, const Box& extent, string filename, bool relax)
// {
//     // read the input file
//     TextInFile in(item, filename, "Tetra vertices");
//     in.addColumn("position x", "length", "pc");
//     in.addColumn("position y", "length", "pc");
//     in.addColumn("position z", "length", "pc");
//     Array coords;
//     while (in.readRow(coords)) _cells.push_back(new Cell(Vec(coords[0], coords[1], coords[2])));
//     in.close();

//     // calculate the Voronoi cells
//     setContext(item);
//     setExtent(extent);
//     buildMesh(relax);
//     buildSearchPerBlock();
// }

////////////////////////////////////////////////////////////////////

// TetraMeshSnapshot::TetraMeshSnapshot(const SimulationItem* item, const Box& extent, SiteListInterface* sli,
//                                          bool relax)
// {
//     // prepare the data
//     int n = sli->numSites();
//     _cells.resize(n);
//     for (int m = 0; m != n; ++m) _cells[m] = new Cell(sli->sitePosition(m));

//     // calculate the Voronoi cells
//     setContext(item);
//     setExtent(extent);
//     buildMesh(relax);
//     buildSearchPerBlock();
// }

////////////////////////////////////////////////////////////////////

TetraMeshSnapshot::TetraMeshSnapshot(const TetraMeshSpatialGrid* grid, const Box& extent)
{
    setContext(grid);
    setExtent(extent);
    buildMesh(grid, true);
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

void TetraMeshSnapshot::buildMesh(const TetraMeshSpatialGrid* grid, bool plc)
{
    tetgenio in, out;
    tetgenbehavior behavior;
    in.firstnumber = 0;

    if (plc)
    {
        in.numberofpoints = 8;
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

        in.numberoffacets = 6;
        // in.numberoffacets += 1;
        in.facetlist = new tetgenio::facet[in.numberoffacets];
        addFacet(&in.facetlist[0], {0, 1, 2, 3});  // Facet 1. bottom
        addFacet(&in.facetlist[1], {4, 5, 6, 7});  // Facet 2. top
        addFacet(&in.facetlist[2], {0, 4, 5, 1});  // Facet 3. front
        addFacet(&in.facetlist[3], {1, 5, 6, 2});  // Facet 4. right
        addFacet(&in.facetlist[4], {2, 6, 7, 3});  // Facet 5. back
        addFacet(&in.facetlist[5], {3, 7, 4, 0});  // Facet 6. left

        behavior.plc = 1;                                           // -p PLC
        behavior.quality = 1;                                       // -q quality mesh
        behavior.fixedvolume = 1;                                   // -a max volume
        behavior.neighout = 2;                                      // -nn neighbors and edges?
        behavior.zeroindex = 1;                                     // -z zero index
        behavior.facesout = 1;                                      // -f faces
        behavior.maxvolume = grid->minVolume() * _extent.volume();  // -a max volume
        
        // behavior.edgesout = 1;     // -e edges
        // behavior.weighted = 1;     // -w weighted
        // behavior.minratio = 5.0;   // -q quality
        // behavior.mindihedral = 5.0;  // -q/ minimal angle

        in.tetunsuitable = [grid](double* pa, double* pb, double* pc, double* pd, double vol) {
            return grid->tetUnsuitable(pa, pb, pc, pd, vol);
        };
    }
    // psc
    else
    {
        in.numberofpoints = 50;
        in.pointlist = new REAL[in.numberofpoints * 3];
        for (int i = 0; i < in.numberofpoints; i++)
        {
            Vec pos = random()->position(_extent);
            in.pointlist[i * 3 + 0] = pos.x();
            in.pointlist[i * 3 + 1] = pos.y();
            in.pointlist[i * 3 + 2] = pos.z();
        }

        behavior.psc = 1;        // -s PSC
        behavior.neighout = 2;   // -nn neighbors and edges?
        behavior.zeroindex = 1;  // -z zero index
        behavior.facesout = 1;   // -f faces
    }

    tetrahedralize(&behavior, &in, &out);

    // tranfser TetGen data in TetraMeshSnapshot data containers
    numTetra = out.numberoftetrahedra;
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
        std::array<Face, 4> faces;

        // calc vertices first
        for (int c = 0; c < 4; c++) vertices[c] = _vertices[out.tetrahedronlist[4 * i + c]];

        for (int c = 0; c < 4; c++)
        {
            // -1 if no neighbor
            int ntetra = out.neighborlist[4 * i + c];

            // find which face is shared with neighbor
            int nface = -1;
            if (ntetra != -1)
            {
                for (int cn = 0; cn < 4; cn++)
                {
                    if (out.neighborlist[4 * ntetra + cn] == i)
                    {
                        nface = cn;
                        break;
                    }
                }
            }

            faces[c] = Face(ntetra, nface);
        }

        _tetrahedra[i] = new Tetra(vertices, faces);

        _centroids.push_back(&_tetrahedra[i]->_centroid);
    }

    // compile statistics
    double minVol = DBL_MAX;
    double maxVol = 0.;
    double totalVol2 = 0.;
    for (int m = 0; m < numTetra; m++)
    {
        double vol = _tetrahedra[m]->volume();
        totalVol2 += vol * vol;
        minVol = min(minVol, vol);
        maxVol = max(maxVol, vol);
    }
    double V = _extent.volume();
    minVol /= V;
    maxVol /= V;
    double avgVol = 1 / (double)numTetra;
    double varVol = (totalVol2 / numTetra / (V * V) - avgVol * avgVol);

    // log neighbor statistics
    log()->info("Done computing tetrahedral tessellation");
    log()->info("  Number of vertices " + std::to_string(numVertices));
    log()->info("  Number of tetrahedra " + std::to_string(numTetra));
    log()->info("  Average volume\% per cell: " + StringUtils::toString(avgVol, 'e'));
    log()->info("  Variance of volume\% per cell: " + StringUtils::toString(varVol, 'e'));
    log()->info("  Minimum volume\% cell: " + StringUtils::toString(minVol, 'e'));
    log()->info("  Maximum volume\% cell: " + StringUtils::toString(maxVol, 'e'));
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
    _rhov.resize(numTetra);
    Array Mv(numTetra);

    // get the maximum temperature, or zero of there is none
    double maxT = useTemperatureCutoff() ? maxTemperature() : 0.;

    // initialize statistics
    double totalOriginalMass = 0;
    double totalMetallicMass = 0;
    double totalEffectiveMass = 0;

    // loop over all sites/cells
    int numIgnored = 0;
    for (int m = 0; m != numTetra; ++m)
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
    if (numTetra) NR::cdf(_cumrhov, Mv);
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
    log()->info("  Number of tetrahedra: " + std::to_string(numTetra));

    // abort if there are no cells
    if (!numTetra) return;

    // construct a single search tree on the site locations of all cells
    log()->info("Building data structure to accelerate searching " + std::to_string(numTetra) + " tetrahedra");
    _blocktrees.resize(1);
    vector<int> ids(numTetra);
    for (int m = 0; m != numTetra; ++m) ids[m] = m;
    _blocktrees[0] = buildTree(ids.begin(), ids.end(), 0);
}

////////////////////////////////////////////////////////////////////

void TetraMeshSnapshot::writeGridPlotFiles(const SimulationItem* probe) const
{
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

    // outputFile.open("data/faces.txt");
    // for (size_t i = 0; i < _tetrahedra.size(); i++)
    // {
    //     const Tetra* tetra = _tetrahedra[i];
    //     bool out = false;
    //     for (int f = 0; f < 4; f++)
    //     {
    //         if (tetra->_faces[f]._ntetra < 0)
    //         {
    //             for (int v : tetra->clockwiseVertices(f))
    //             {
    //                 const Vec* r = tetra->_vertices[v];
    //                 outputFile << r->x() << ", " << r->y() << ", " << r->z() << "\n";
    //             }
    //         }
    //     }
    // }
    // outputFile.close();

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

    // int test = 0;

    // look for the closest centroid in this block using the search tree
    Node* tree = _blocktrees[b];
    if (tree)
    {
        // /*
        // full traversal algorithm
        int enteringFace = -1;
        int leavingFace = -1;
        int m = tree->m();

        const Tetra* tetra = _tetrahedra[m];

        Vec pos = tetra->_centroid;
        Direction dir(bfr - tetra->_centroid);
        double dist = dir.norm();  // keep subtracting ds until dist < 0
        dir /= dist;

        while (true)
        {
            tetra = _tetrahedra[m];

            // find entering face using a single Plücker product
            if (enteringFace == -1)
            {
                enteringFace = findEnteringFace(tetra, pos, dir);
                if (enteringFace == -1) break;
            }

            // the translated Plücker moment in the local coordinate system
            Vec moment = Vec::cross(dir, pos - *tetra->_vertices[enteringFace]);
            std::array<int, 3> cv = Tetra::clockwiseVertices(enteringFace);

            // 2 step decision tree
            double prod0 = Vec::dot(moment, tetra->getEdge(cv[0], enteringFace));
            int clock0 = prod0 < 0;
            // if clockwise move clockwise else move cclockwise
            int i = clock0 ? 1 : 2;
            double prodi = Vec::dot(moment, tetra->getEdge(cv[i], enteringFace));
            int cclocki = prodi >= 0;

            double ds = DBL_MAX;

            // use plane intersection algorithm if Plücker products are ambiguous
            if (prod0 == 0. || prodi == 0.)
            {
                for (int face : cv)
                {
                    const Vec& n = tetra->_faces[face]._normal;
                    double ndotk = Vec::dot(n, dir);
                    if (ndotk > 0)
                    {
                        const Vec& v = *tetra->_vertices[enteringFace];
                        double dq = Vec::dot(n, v - pos) / ndotk;
                        if (dq < ds)
                        {
                            ds = dq;
                            leavingFace = face;
                        }
                    }
                }
            }
            // use Maria (2017) algorithm otherwise
            else
            {
                // decision table for clock0 and cclocki
                // 1 1 -> 2
                // 0 0 -> 1
                // 1 0 -> 0
                // 0 1 -> 0
                constexpr int dtable[2][2] = {{1, 0}, {0, 2}};
                leavingFace = cv[dtable[clock0][cclocki]];
                const Vec& n = tetra->_faces[leavingFace]._normal;
                const Vec& v = *tetra->_vertices[enteringFace];
                double ndotk = Vec::dot(n, dir);

                ds = Vec::dot(n, v - pos) / ndotk;
            }

            // if ds is too close to leaving face we recalculate cellIndex to avoid traversing when ds ~ 0
            // this might actually slow down the traversal

            // if no exit point was found, advance the current point by a small distance,
            // recalculate the cell index, and return to the start of the loop
            if (leavingFace == -1)
            {
                break;  // traversal failed
            }
            // otherwise set the current point to the exit point and return the path segment
            else
            {
                pos += ds * dir;
                dist -= ds;
                if (dist <= 0) return m;

                m = tetra->_faces[leavingFace]._ntetra;
                enteringFace = tetra->_faces[leavingFace]._nface;

                if (m < 0) break;
            }
        }
        // */

        /*
        // traverse towards input point from nearest centroid
        int m = tree->nearest(bfr, _centroids)->m();
        const Tetra* root = _tetrahedra[m];

        test++;
        if (root->inside(bfr)) return test;

        Vec pos = root->_centroid;
        Vec dir = bfr - pos;
        double dist = dir.norm();
        dir /= dist;

        // start by checking all faces
        std::vector<int> faces = {0, 1, 2, 3};

        while (true)
        {
            const Tetra* tetra = _tetrahedra[m];

            double ds = DBL_MAX;
            int leavingFace = -1;

            // find exit face
            for (int face : faces)
            {
                const Vec& n = tetra->_faces[face]._normal;
                double ndotk = Vec::dot(n, dir);
                if (ndotk > 0)
                {
                    const Vec& v = *tetra->_vertices[(face + 1) % 4];
                    double dq = Vec::dot(n, v - pos) / ndotk;
                    if (dq < ds)
                    {
                        ds = dq;
                        leavingFace = face;
                    }
                }
            }

            if (leavingFace == -1) break;

            // propagate and check inside if we passed the input point
            pos += ds * dir;
            dist -= ds;
            if (dist <= 0)
            {
                test++;
                if (tetra->inside(bfr)) return test;
            }

            m = tetra->_faces[leavingFace]._ntetra;
            int enteringFace = tetra->_faces[leavingFace]._nface;
            faces = {(enteringFace + 1) % 4, (enteringFace + 2) % 4, (enteringFace + 3) % 4};

            if (m < 0) break;
        }
        */

        /*
        // use a breadth-first search over the neighboring tetrahedra
        int root = tree->nearest(bfr, _centroids)->m();
        std::queue<int> queue;
        std::unordered_set<int> explored;
        queue.push(root);

        // Vec dir = bfr - _tetrahedra[root]->_centroid;

        int t;
        while (queue.size() > 0)
        {
            t = queue.front();
            queue.pop();
            const Tetra* tetra = _tetrahedra[t];

            test++;
            if (tetra->inside(bfr)) return test;
            explored.insert(t);

            Vec dir = bfr - tetra->_centroid;
            for (const Face& face : tetra->_faces)
            {
                int n = face._ntetra;
                if (n != -1 && (explored.find(n) == explored.end()))  // if not already explored
                {
                    if (Vec::dot(dir, face._normal) > 0)
                    {
                        queue.push(n);
                        explored.insert(n);
                    }
                }
            }
        }
        */

        log()->error("search tree failed to find the tetrahedron");  // change to warning?
    }

    // if there is no search tree or search tree failed, simply loop over all tetrahedra in the block
    for (int t : _blocklists[b])
    {
        if (_tetrahedra[t]->inside(bfr)) return t;
    }

    log()->error("cellIndex failed to find the tetrahedron");  // change to warning?
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
    int _mr;
    int _enteringFace;

public:
    MySegmentGenerator(const TetraMeshSnapshot* grid) : _grid(grid) {}

    // #define WRITE

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
            _enteringFace = -1;

            // if the photon packet started outside the grid, return the corresponding nonzero-length segment;
            // otherwise fall through to determine the first actual segment
            if (ds() > 0.) return true;

#ifdef WRITE
            startwriting();
            write(0, _mr);
#endif
        }

        // intentionally falls through
        if (state() == State::Inside)
        {
            int leavingFace = -1;
            // loop in case no exit point was found (which should happen only rarely)
            while (true)
            {
                const Tetra* tetra = _grid->_tetrahedra[_mr];
                Position pos = r();
                Direction dir = k();

                // find entering face using a single Plücker product
                if (_enteringFace == -1)
                {
                    _enteringFace = findEnteringFace(tetra, pos, dir);
                }

                // the translated Plücker moment in the local coordinate system
                Vec moment = Vec::cross(dir, pos - *tetra->_vertices[_enteringFace]);
                std::array<int, 3> cv = Tetra::clockwiseVertices(_enteringFace);

                // 2 step decision tree
                double prod0 = Vec::dot(moment, tetra->getEdge(cv[0], _enteringFace));
                int clock0 = prod0 < 0;
                // if clockwise move clockwise else move cclockwise
                int i = clock0 ? 1 : 2;
                double prodi = Vec::dot(moment, tetra->getEdge(cv[i], _enteringFace));
                int cclocki = prodi >= 0;

                double ds = DBL_MAX;

                // use plane intersection algorithm if Plücker products are ambiguous
                if (prod0 == 0. || prodi == 0.)
                {
                    for (int face : cv)
                    {
                        const Vec& n = tetra->_faces[face]._normal;
                        double ndotk = Vec::dot(n, dir);
                        if (ndotk > 0)
                        {
                            const Vec& v = *tetra->_vertices[_enteringFace];
                            double dq = Vec::dot(n, v - pos) / ndotk;
                            if (dq < ds)
                            {
                                ds = dq;
                                leavingFace = face;
                            }
                        }
                    }
                }
                // use Maria (2017) algorithm otherwise
                else
                {
                    // decision table for clock0 and cclocki
                    // 1 1 -> 2
                    // 0 0 -> 1
                    // 1 0 -> 0
                    // 0 1 -> 0
                    constexpr int dtable[2][2] = {{1, 0}, {0, 2}};
                    leavingFace = cv[dtable[clock0][cclocki]];
                    const Vec& n = tetra->_faces[leavingFace]._normal;
                    const Vec& v = *tetra->_vertices[_enteringFace];
                    double ndotk = Vec::dot(n, dir);

                    // if ndotk == 0. we are inside the leaving face, we should just move to the neighboring tetra and not move (ds=0)
                    // once it moved to its neighbor, it can't move back to this cell and it will traverse the next step correctly
                    // if (ndotk == 0.)
                    //     ds = 0;
                    // else
                    //     ds = Vec::dot(n, v - pos) / ndotk;

                    ds = Vec::dot(n, v - pos) / ndotk;
                }

                // if ds is too close to leaving face we recalculate cellIndex to avoid traversing when ds ~ 0
                // this might actually slow down the traversal

                // if no exit point was found, advance the current point by a small distance,
                // recalculate the cell index, and return to the start of the loop
                if (leavingFace == -1 || ds < _grid->_eps)
                {
                    // move towards centroid?
                    propagater(_grid->_eps);
                    _mr = _grid->cellIndex(r());

                    if (_mr < 0)
                    {
                        setState(State::Outside);
                        return false;
                    }

                    _enteringFace = -1;
#ifdef WRITE
                    write(ds, -1);
#endif
                }
                // otherwise set the current point to the exit point and return the path segment
                else
                {
                    propagater(ds);
                    setSegment(_mr, ds);
                    _mr = tetra->_faces[leavingFace]._ntetra;

                    if (_mr < 0)
                    {
                        setState(State::Outside);
                        return false;
                    }

                    _enteringFace = tetra->_faces[leavingFace]._nface;

#ifdef WRITE
                    write(ds, leavingFace);
#endif

                    return true;
                }
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
