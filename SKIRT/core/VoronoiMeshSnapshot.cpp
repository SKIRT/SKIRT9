/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "VoronoiMeshSnapshot.hpp"
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

// class to hold the information about a Voronoi cell that is relevant for calculating paths and densities
class VoronoiMeshSnapshot::Cell : public Box  // enclosing box
{
private:
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

    // initializes the receiver with information taken from the specified fully computed Voronoi cell
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

    // clears the information taken from a Voronoi cell so it can be reinitialized
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

    // writes the Voronoi cell geometry to the serialized data buffer, preceded by the specified cell index,
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

    // reads the Voronoi cell geometry from the serialized data buffer
    void readGeometry(SerializedRead& rdata)
    {
        rdata.read(*this);  // extent
        rdata.read(_c);
        rdata.read(_volume);
        rdata.read(_neighbors);
    }
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
}

////////////////////////////////////////////////////////////////////

// class to hold a node in the binary search tree (see en.wikipedia.org/wiki/Kd-tree)
class VoronoiMeshSnapshot::Node
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

VoronoiMeshSnapshot::VoronoiMeshSnapshot() {}

////////////////////////////////////////////////////////////////////

VoronoiMeshSnapshot::~VoronoiMeshSnapshot()
{
    for (auto cell : _cells) delete cell;
    for (auto tree : _blocktrees) delete tree;
}

////////////////////////////////////////////////////////////////////

void VoronoiMeshSnapshot::readAndClose()
{
    // read the site info into memory
    Array prop;
    while (infile()->readRow(prop)) _cells.push_back(new Cell(prop));

    // close the file
    Snapshot::readAndClose();

    // if we are allowed to build a Voronoi mesh
    if (!_foregoVoronoiMesh)
    {
        // calculate the Voronoi cells
        buildMesh(false);

        // if a mass density policy has been set, calculate masses and densities and build the search data structure
        if (hasMassDensityPolicy()) calculateDensityAndMass();
        if (hasMassDensityPolicy() || needGetEntities()) buildSearchPerBlock();
    }

    // if we forego building a Voronoi mesh, there is a density policy by definition
    else
    {
        calculateVolume();
        calculateDensityAndMass();
        buildSearchSingle();
    }
}

////////////////////////////////////////////////////////////////////

void VoronoiMeshSnapshot::setExtent(const Box& extent)
{
    _extent = extent;
    _eps = 1e-12 * extent.widths().norm();
}

////////////////////////////////////////////////////////////////////

void VoronoiMeshSnapshot::foregoVoronoiMesh()
{
    _foregoVoronoiMesh = true;
}

////////////////////////////////////////////////////////////////////

VoronoiMeshSnapshot::VoronoiMeshSnapshot(const SimulationItem* item, const Box& extent, string filename, bool relax)
{
    // read the input file
    TextInFile in(item, filename, "Voronoi sites");
    in.addColumn("position x", "length", "pc");
    in.addColumn("position y", "length", "pc");
    in.addColumn("position z", "length", "pc");
    Array coords;
    while (in.readRow(coords)) _cells.push_back(new Cell(Vec(coords[0], coords[1], coords[2])));
    in.close();

    // calculate the Voronoi cells
    setContext(item);
    setExtent(extent);
    buildMesh(relax);
    buildSearchPerBlock();
}

////////////////////////////////////////////////////////////////////

VoronoiMeshSnapshot::VoronoiMeshSnapshot(const SimulationItem* item, const Box& extent, SiteListInterface* sli,
                                         bool relax)
{
    // prepare the data
    int n = sli->numSites();
    _cells.resize(n);
    for (int m = 0; m != n; ++m) _cells[m] = new Cell(sli->sitePosition(m));

    // calculate the Voronoi cells
    setContext(item);
    setExtent(extent);
    buildMesh(relax);
    buildSearchPerBlock();
}

////////////////////////////////////////////////////////////////////

VoronoiMeshSnapshot::VoronoiMeshSnapshot(const SimulationItem* item, const Box& extent, const vector<Vec>& sites,
                                         bool relax)
{
    // prepare the data
    int n = sites.size();
    _cells.resize(n);
    for (int m = 0; m != n; ++m) _cells[m] = new Cell(sites[m]);

    // calculate the Voronoi cells
    setContext(item);
    setExtent(extent);
    buildMesh(relax);
    buildSearchPerBlock();
}

////////////////////////////////////////////////////////////////////

namespace
{
    // maximum number of Voronoi sites processed between two invocations of infoIfElapsed()
    const int logProgressChunkSize = 1000;

    // maximum number of Voronoi grid construction iterations
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

void VoronoiMeshSnapshot::buildMesh(bool relax)
{
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

    // repeat grid construction until none of the cells have zero volume
    int numIterations = 0;
    while (true)
    {
        // add the final sites to a temporary Voronoi container, using the cell index m as ID
        voro::container vcon(_extent.xmin(), _extent.xmax(), _extent.ymin(), _extent.ymax(), _extent.zmin(),
                             _extent.zmax(), _nb, _nb, _nb, false, false, false, 16);
        for (int m = 0; m != numCells; ++m)
        {
            Vec r = _cells[m]->position();
            vcon.put(m, r.x(), r.y(), r.z());
        }

        // for each site:
        //   - compute the corresponding cell in the Voronoi tesselation
        //   - extract and copy the relevant information to the cell object with the corresponding index in our vector
        log()->info("Constructing Voronoi tessellation with " + std::to_string(numCells) + " cells");
        log()->infoSetElapsed(numCells);
        auto parallel = log()->find<ParallelFactory>()->parallelDistributed();
        parallel->call(numCells, [this, &vcon](size_t firstIndex, size_t numIndices) {
            // allocate a separate cell calculator for each thread to avoid conflicts
            voro::voro_compute<voro::container> vcompute(vcon, _nb, _nb, _nb);
            // allocate space for the resulting cell info
            voro::voronoicell_neighbor vcell;

            // loop over all cells and work on the ones that have a particle index in our dedicated range
            // (we cannot access cells in the container based on cell index m without building an extra data structure)
            int numDone = 0;
            voro::c_loop_all vloop(vcon);
            if (vloop.start()) do
                {
                    size_t m = vloop.pid();
                    if (m >= firstIndex && m < firstIndex + numIndices)
                    {
                        // compute the cell and copy all relevant information to the cell object that will stay around
                        bool ok = vcompute.compute_cell(vcell, vloop.ijk, vloop.q, vloop.i, vloop.j, vloop.k);
                        if (ok) _cells[m]->init(vcell);

                        // log message if the minimum time has elapsed
                        numDone = (numDone + 1) % logProgressChunkSize;
                        if (numDone == 0) log()->infoIfElapsed("Computed Voronoi cells: ", logProgressChunkSize);
                    }
                } while (vloop.inc());
            if (numDone > 0) log()->infoIfElapsed("Computed Voronoi cells: ", numDone);
        });

        // communicate the calculated cell information between parallel processes, if needed
        if (ProcessManager::isMultiProc())
        {
            auto producer = [this](vector<double>& data) {
                SerializedWrite wdata(data);
                int numCells = _cells.size();
                for (int m = 0; m != numCells; ++m) _cells[m]->writeGeometryIfPresent(wdata, m);
            };
            auto consumer = [this](const vector<double>& data) {
                SerializedRead rdata(data);
                while (!rdata.empty()) _cells[rdata.readInt()]->readGeometry(rdata);
            };
            ProcessManager::broadcastAllToAll(producer, consumer);
        }

        // discover invalid cells with zero volume and/or with neighbors that are not mutual
        log()->info("Verifying Voronoi tessellation");
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
            throw FATALERROR("Still " + std::to_string(invalid.size()) + " invalid Voronoi cells after "
                             + std::to_string(maxConstructionIterations)
                             + " iterations of constructing the tessellation");
        }

        // remove invalid cells and prepare to repeat construction
        log()->warning("Removing sites for " + std::to_string(invalid.size())
                       + " invalid Voronoi cells and reconstructing the tessellation");
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
    log()->info("Done computing Voronoi tessellation with " + std::to_string(numCells) + " cells");
    log()->info("  Average number of neighbors per cell: " + StringUtils::toString(avgNeighbors, 'f', 1));
    log()->info("  Minimum number of neighbors per cell: " + std::to_string(minNeighbors));
    log()->info("  Maximum number of neighbors per cell: " + std::to_string(maxNeighbors));
}

////////////////////////////////////////////////////////////////////

void VoronoiMeshSnapshot::calculateVolume()
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

void VoronoiMeshSnapshot::calculateDensityAndMass()
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

VoronoiMeshSnapshot::Node* VoronoiMeshSnapshot::buildTree(vector<int>::iterator first, vector<int>::iterator last,
                                                          int depth) const
{
    auto length = last - first;
    if (length > 0)
    {
        auto median = length >> 1;
        std::nth_element(first, first + median, last, [this, depth](int m1, int m2) {
            return m1 != m2 && lessthan(_cells[m1]->position(), _cells[m2]->position(), depth % 3);
        });
        return new VoronoiMeshSnapshot::Node(*(first + median), depth, buildTree(first, first + median, depth + 1),
                                             buildTree(first + median + 1, last, depth + 1));
    }
    return nullptr;
}

////////////////////////////////////////////////////////////////////

void VoronoiMeshSnapshot::buildSearchPerBlock()
{
    // abort if there are no cells
    int numCells = _cells.size();
    if (!numCells) return;

    log()->info("Building data structures to accelerate searching the Voronoi tesselation");

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

void VoronoiMeshSnapshot::buildSearchSingle()
{
    // log the number of sites
    int numCells = _cells.size();
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

bool VoronoiMeshSnapshot::isPointClosestTo(Vec r, int m, const vector<int>& ids) const
{
    double target = _cells[m]->squaredDistanceTo(r);
    for (int id : ids)
    {
        if (id >= 0 && _cells[id]->squaredDistanceTo(r) < target) return false;
    }
    return true;
}

////////////////////////////////////////////////////////////////////

void VoronoiMeshSnapshot::writeGridPlotFiles(const SimulationItem* probe) const
{
    // create the plot files
    SpatialGridPlotFile plotxy(probe, probe->itemName() + "_grid_xy");
    SpatialGridPlotFile plotxz(probe, probe->itemName() + "_grid_xz");
    SpatialGridPlotFile plotyz(probe, probe->itemName() + "_grid_yz");
    SpatialGridPlotFile plotxyz(probe, probe->itemName() + "_grid_xyz");

    // load all sites in a Voro container
    int numCells = _cells.size();
    voro::container vcon(_extent.xmin(), _extent.xmax(), _extent.ymin(), _extent.ymax(), _extent.zmin(), _extent.zmax(),
                         _nb, _nb, _nb, false, false, false, 16);
    for (int m = 0; m != numCells; ++m)
    {
        Vec r = _cells[m]->position();
        vcon.put(m, r.x(), r.y(), r.z());
    }

    // for each site, compute the corresponding cell and output its edges
    log()->info("Writing plot files for Voronoi tessellation with " + std::to_string(numCells) + " cells");
    log()->infoSetElapsed(numCells);
    voro::voronoicell_neighbor vcell;
    int numDone = 0;
    voro::c_loop_all vloop(vcon);
    if (vloop.start()) do
        {
            // compute the cell
            vcon.compute_cell(vcell, vloop);

            // get the edges of the cell
            double x, y, z;
            vloop.pos(x, y, z);
            vector<double> coords;
            vcell.vertices(x, y, z, coords);
            vector<int> indices;
            vcell.face_vertices(indices);

            // write the edges of the cell to the plot files
            Box bounds = _cells[vloop.pid()]->extent();
            if (bounds.zmin() <= 0 && bounds.zmax() >= 0) plotxy.writePolyhedron(coords, indices);
            if (bounds.ymin() <= 0 && bounds.ymax() >= 0) plotxz.writePolyhedron(coords, indices);
            if (bounds.xmin() <= 0 && bounds.xmax() >= 0) plotyz.writePolyhedron(coords, indices);
            if (vloop.pid() <= 1000) plotxyz.writePolyhedron(coords, indices);

            // log message if the minimum time has elapsed
            numDone++;
            if (numDone % 2000 == 0) log()->infoIfElapsed("Computed Voronoi cells: ", 2000);
        } while (vloop.inc());
}

////////////////////////////////////////////////////////////////////

Box VoronoiMeshSnapshot::extent() const
{
    return _extent;
}

////////////////////////////////////////////////////////////////////

int VoronoiMeshSnapshot::numEntities() const
{
    return _cells.size();
}

////////////////////////////////////////////////////////////////////

Position VoronoiMeshSnapshot::position(int m) const
{
    return Position(_cells[m]->position());
}

////////////////////////////////////////////////////////////////////

Position VoronoiMeshSnapshot::centroidPosition(int m) const
{
    return Position(_cells[m]->centroid());
}

////////////////////////////////////////////////////////////////////

double VoronoiMeshSnapshot::volume(int m) const
{
    return _cells[m]->volume();
}

////////////////////////////////////////////////////////////////////

Box VoronoiMeshSnapshot::extent(int m) const
{
    return _cells[m]->extent();
}

////////////////////////////////////////////////////////////////////

double VoronoiMeshSnapshot::density(int m) const
{
    return _rhov[m];
}

////////////////////////////////////////////////////////////////////

double VoronoiMeshSnapshot::density(Position bfr) const
{
    int m = cellIndex(bfr);
    return m >= 0 ? _rhov[m] : 0;
}

////////////////////////////////////////////////////////////////////

double VoronoiMeshSnapshot::mass() const
{
    return _mass;
}

////////////////////////////////////////////////////////////////////

Position VoronoiMeshSnapshot::generatePosition(int m) const
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

Position VoronoiMeshSnapshot::generatePosition() const
{
    // if there are no sites, return the origin
    if (_cells.empty()) return Position();

    // select a site according to its mass contribution
    int m = NR::locateClip(_cumrhov, random()->uniform());

    return generatePosition(m);
}

////////////////////////////////////////////////////////////////////

int VoronoiMeshSnapshot::cellIndex(Position bfr) const
{
    // make sure the position is inside the domain
    if (!_extent.contains(bfr)) return -1;

    // determine the block in which the point falls
    // if we didn't build a Voronoi mesh, the search tree is always in the first "block"
    int b = 0;
    if (!_foregoVoronoiMesh)
    {
        int i, j, k;
        _extent.cellIndices(i, j, k, bfr, _nb, _nb, _nb);
        b = i * _nb2 + j * _nb + k;
    }

    // look for the closest site in this block, using the search tree if there is one
    Node* tree = _blocktrees[b];
    if (tree) return tree->nearest(bfr, _cells)->m();

    // if there is no search tree, simply loop over the index list
    const vector<int>& ids = _blocklists[b];
    int m = -1;
    double mdist = DBL_MAX;
    int n = ids.size();
    for (int i = 0; i < n; i++)
    {
        double idist = _cells[ids[i]]->squaredDistanceTo(bfr);
        if (idist < mdist)
        {
            m = ids[i];
            mdist = idist;
        }
    }
    return m;
}

////////////////////////////////////////////////////////////////////

const Array& VoronoiMeshSnapshot::properties(int m) const
{
    return _cells[m]->properties();
}

////////////////////////////////////////////////////////////////////

int VoronoiMeshSnapshot::nearestEntity(Position bfr) const
{
    return _blocktrees.size() ? cellIndex(bfr) : -1;
}

////////////////////////////////////////////////////////////////////

void VoronoiMeshSnapshot::getEntities(EntityCollection& entities, Position bfr) const
{
    entities.addSingle(cellIndex(bfr));
}

////////////////////////////////////////////////////////////////////

class VoronoiMeshSnapshot::MySegmentGenerator : public PathSegmentGenerator
{
    const VoronoiMeshSnapshot* _grid{nullptr};
    int _mr{-1};

public:
    MySegmentGenerator(const VoronoiMeshSnapshot* grid) : _grid(grid) {}

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

void VoronoiMeshSnapshot::getEntities(EntityCollection& entities, Position bfr, Direction bfk) const
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

std::unique_ptr<PathSegmentGenerator> VoronoiMeshSnapshot::createPathSegmentGenerator() const
{
    return std::make_unique<MySegmentGenerator>(this);
}

////////////////////////////////////////////////////////////////////
