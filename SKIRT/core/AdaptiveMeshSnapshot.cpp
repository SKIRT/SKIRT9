/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "AdaptiveMeshSnapshot.hpp"
#include "EntityCollection.hpp"
#include "FatalError.hpp"
#include "Log.hpp"
#include "NR.hpp"
#include "PathSegmentGenerator.hpp"
#include "Random.hpp"
#include "SpatialGridPath.hpp"
#include "StringUtils.hpp"
#include "TextInFile.hpp"
#include "Units.hpp"

////////////////////////////////////////////////////////////////////

/* The Node class is a helper class used to represent individual nodes in a tree data
   structure. An Node instance can represent a leaf or a nonleaf node. A nonleaf node
   maintains a list of pointers to its children. A leaf node instead keeps a list of pointers
   to the most likely neighbor for each of the six cell walls. */
class AdaptiveMeshSnapshot::Node : public Box
{
public:
    /* The constructor receives the node's extent, and reads the other node data from the
       following line in the specified input file. It then recursively constructs any child
       nodes. In addition to constructing the new node(s), the constructor also adds leaf node
       pointers to the vector held by the AdaptiveMeshSnapshot class. */
    Node(const Box& extent, TextInFile* infile, vector<Node*>& leafnodes) : Box(extent)
    {
        // if this is a nonleaf line, process it
        if (infile->readNonLeaf(_Nx, _Ny, _Nz))
        {
            _m = -1;

            // construct and store our children, in local Morton order
            _nodes.resize(_Nx * _Ny * _Nz);
            int m = 0;
            for (int k = 0; k < _Nz; k++)
                for (int j = 0; j < _Ny; j++)
                    for (int i = 0; i < _Nx; i++)
                    {
                        Vec r0 = extent.fracPos(i, j, k, _Nx, _Ny, _Nz);
                        Vec r1 = extent.fracPos(i + 1, j + 1, k + 1, _Nx, _Ny, _Nz);
                        _nodes[m++] = new Node(Box(r0, r1), infile, leafnodes);
                    }
        }

        // if this is not a nonleaf line, it should be a leaf line
        else
        {
            _Nx = 0;
            _Ny = 0;
            _Nz = 0;

            // read a leaf line and detect premature end-of file
            if (!infile->readRow(_properties))
                throw FATALERROR("Reached end of file in adaptive mesh data before all nodes were read");

            // add this leaf node to the list
            _m = leafnodes.size();
            leafnodes.push_back(this);
        }
    }

    /* This function adds neighbor information to the receiving leaf node. Specifically, it
       constructs a list of the node's most likely neighbor at each of its six walls. The function
       does nothing if neighbor information has already been added, or if the node is a nonleaf
       node. The first argument specifies the adaptive mesh root node, which is queried by this
       function to locate the node's neighbors. The second argument specifies a very small offset
       (relative to the domain size) used to determine a location just beyond a node wall. */
    void addNeighbors(Node* root, double eps)
    {
        // skip if this node has children or if neighbors already have been added
        if (_nodes.empty())
        {
            // node center
            Vec ctr = center();

            // determine the node just beyond the center of each wall (or null pointer for domain walls)
            _nodes.resize(6);
            _nodes[BACK] = root->leaf(ctr + Vec(-eps, 0, 0));
            _nodes[FRONT] = root->leaf(ctr + Vec(+eps, 0, 0));
            _nodes[LEFT] = root->leaf(ctr + Vec(0, -eps, 0));
            _nodes[RIGHT] = root->leaf(ctr + Vec(0, +eps, 0));
            _nodes[BOTTOM] = root->leaf(ctr + Vec(0, 0, -eps));
            _nodes[TOP] = root->leaf(ctr + Vec(0, 0, +eps));
        }
    }

    /* The destructor releases the node's children if it is a nonleaf node. */
    ~Node()
    {
        if (!isLeaf())
            for (auto node : _nodes) delete node;
    }

    /* For leaf nodes, this function returns the Morton order cell index of the corresponding
       cell. For nonleaf nodes, this function returns -1. */
    int cellIndex() const { return _m; }

    /* This function returns true if the node is a leaf node, false if it is a nonleaf node. */
    bool isLeaf() const { return _m >= 0; }

    /* This function returns a pointer to the node's immediate child that contains the specified
       point, assuming that the point is inside the node (which is not verified). This function
       causes undefined behavior if the node is a leaf node. */
    const Node* child(Vec r) const
    {
        // estimate the child node indices; this may be off by one due to rounding errors
        int i, j, k;
        cellIndices(i, j, k, r, _Nx, _Ny, _Nz);

        // get the estimated node using local Morton order
        const Node* node = _nodes[(k * _Ny + j) * _Nx + i];

        // if the point is NOT in the node, correct the indices and get the new node
        if (!node->contains(r))
        {
            if (r.x() < node->xmin())
                i--;
            else if (r.x() > node->xmax())
                i++;
            if (r.y() < node->ymin())
                j--;
            else if (r.y() > node->ymax())
                j++;
            if (r.z() < node->zmin())
                k--;
            else if (r.z() > node->zmax())
                k++;
            node = _nodes[(k * _Ny + j) * _Nx + i];
            if (!node->contains(r)) throw FATALERROR("Can't locate the appropriate child node");
        }
        return node;
    }

    /* This function returns a pointer to the deepest node in the child hierarchy of this node
       that contains the specified point, or null if the point is outside the node. It uses the
       child() function repeatedly to locate the appropriate node. */
    const Node* leaf(Vec r) const
    {
        if (!contains(r)) return nullptr;

        const Node* node = this;
        while (!node->isLeaf()) node = node->child(r);
        return node;
    }

    /* This enum lists a constant for each of the walls in a node. The x-coordinate increases
       from BACK to FRONT, the y-coordinate increases from LEFT to RIGHT, and the z-coordinate
       increases from BOTTOM to TOP. */
    enum Wall { BACK = 0, FRONT, LEFT, RIGHT, BOTTOM, TOP };

    /* This function returns a pointer to the leaf node just beyond a given wall that contains
       the specified position, or null if the specified position is not inside the most likely
       neighbor for that wall, or if neighbors have not been added for the node, or if this is
       not a leaf node. */
    const Node* neighbor(Wall wall, Vec r) const
    {
        const Node* node = _nodes.size() == 6 ? _nodes[wall] : nullptr;
        if (node && node->contains(r))
            return node;
        else
            return nullptr;
    }

    /* This function returns the user properties for leaf nodes; for nonleaf nodes the
       returned array is empty. */
    const Array& properties() { return _properties; }

private:
    int _Nx, _Ny, _Nz;           // number of grid cells in each direction; zero for leaf nodes
    int _m;                      // Morton order index for the cell represented by this leaf node; -1 for nonleaf nodes
    vector<const Node*> _nodes;  // pointers to children (nonleaf nodes) or neighbors (leaf nodes)
    Array _properties;           // user-defined properties
};

////////////////////////////////////////////////////////////////////

AdaptiveMeshSnapshot::AdaptiveMeshSnapshot() {}

////////////////////////////////////////////////////////////////////

AdaptiveMeshSnapshot::~AdaptiveMeshSnapshot()
{
    delete _root;
}

////////////////////////////////////////////////////////////////////

void AdaptiveMeshSnapshot::readAndClose()
{
    // construct the root node, and recursively all other nodes;
    // this also fills the _cells vector
    _root = new Node(_extent, infile(), _cells);

    // verify that all data was read and close the file
    Array dummy;
    if (infile()->readRow(dummy)) throw FATALERROR("Superfluous lines in adaptive mesh data after all nodes were read");
    Snapshot::readAndClose();

    // log nr of cells
    log()->info("  Number of leaf cells: " + std::to_string(_cells.size()));

    // if a mass density policy has been set, calculate masses and densities for all cells
    if (hasMassDensityPolicy())
    {
        // allocate vectors for mass and density
        size_t n = _cells.size();
        Array Mv(n);
        _rhov.resize(n);

        // get the maximum temperature, or zero of there is none
        double maxT = useTemperatureCutoff() ? maxTemperature() : 0.;

        // initialize statistics
        double totalOriginalMass = 0;
        double totalMetallicMass = 0;
        double totalEffectiveMass = 0;

        // loop over all leaf cells
        int numIgnored = 0;
        for (size_t m = 0; m != n; ++m)
        {
            const Array& prop = _cells[m]->properties();

            // original mass is zero if temperature is above cutoff or if imported mass/density is not positive
            double originalMass = 0.;
            if (maxT && prop[temperatureIndex()] > maxT)
                numIgnored++;
            else
                originalMass =
                    max(0., massIndex() >= 0 ? prop[massIndex()] : prop[densityIndex()] * _cells[m]->volume());

            double metallicMass = originalMass * (useMetallicity() ? prop[metallicityIndex()] : 1.);
            double effectiveMass = metallicMass * multiplier();

            Mv[m] = effectiveMass;
            _rhov[m] = effectiveMass / _cells[m]->volume();

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
}

////////////////////////////////////////////////////////////////////

void AdaptiveMeshSnapshot::setExtent(const Box& extent)
{
    _extent = extent;
    _eps = 1e-12 * extent.widths().norm();
}

////////////////////////////////////////////////////////////////////

void AdaptiveMeshSnapshot::addNeighbors()
{
    for (auto cell : _cells) cell->addNeighbors(_root, _eps);
}

////////////////////////////////////////////////////////////////////

Box AdaptiveMeshSnapshot::extent() const
{
    return _extent;
}

////////////////////////////////////////////////////////////////////

int AdaptiveMeshSnapshot::numEntities() const
{
    return _cells.size();
}

////////////////////////////////////////////////////////////////////

Position AdaptiveMeshSnapshot::position(int m) const
{
    return Position(_cells[m]->center());
}

////////////////////////////////////////////////////////////////////

double AdaptiveMeshSnapshot::volume(int m) const
{
    return _cells[m]->volume();
}

////////////////////////////////////////////////////////////////////

double AdaptiveMeshSnapshot::diagonal(int m) const
{
    return _cells[m]->diagonal();
}

////////////////////////////////////////////////////////////////////

Box AdaptiveMeshSnapshot::extent(int m) const
{
    return _cells[m]->extent();
}

////////////////////////////////////////////////////////////////////

double AdaptiveMeshSnapshot::density(int m) const
{
    return _rhov[m];
}

////////////////////////////////////////////////////////////////////

double AdaptiveMeshSnapshot::density(Position bfr) const
{
    int m = cellIndex(bfr);
    return m >= 0 ? _rhov[m] : 0;
}

////////////////////////////////////////////////////////////////////

double AdaptiveMeshSnapshot::mass() const
{
    return _mass;
}

////////////////////////////////////////////////////////////////////

Position AdaptiveMeshSnapshot::generatePosition(int m) const
{
    return random()->position(_cells[m]->extent());
}

////////////////////////////////////////////////////////////////////

Position AdaptiveMeshSnapshot::generatePosition() const
{
    // if there are no cells, return the origin
    if (_cells.empty()) return Position();

    // select a cell according to its mass contribution
    int m = NR::locateClip(_cumrhov, random()->uniform());

    return generatePosition(m);
}

////////////////////////////////////////////////////////////////////

int AdaptiveMeshSnapshot::cellIndex(Position bfr) const
{
    const Node* cell = _root->leaf(bfr);
    return cell ? cell->cellIndex() : -1;
}

////////////////////////////////////////////////////////////////////

const Array& AdaptiveMeshSnapshot::properties(int m) const
{
    return _cells[m]->properties();
}

////////////////////////////////////////////////////////////////////

int AdaptiveMeshSnapshot::nearestEntity(Position bfr) const
{
    return cellIndex(bfr);
}

////////////////////////////////////////////////////////////////////

void AdaptiveMeshSnapshot::getEntities(EntityCollection& entities, Position bfr) const
{
    entities.addSingle(cellIndex(bfr));
}

////////////////////////////////////////////////////////////////////

class AdaptiveMeshSnapshot::MySegmentGenerator : public PathSegmentGenerator
{
    const AdaptiveMeshSnapshot* _grid{nullptr};
    const Node* _node{nullptr};

public:
    MySegmentGenerator(const AdaptiveMeshSnapshot* grid) : _grid(grid) {}

    bool next() override
    {
        switch (state())
        {
            case State::Unknown:
            {
                // try moving the photon packet inside the grid; if this is impossible, return an empty path
                if (!moveInside(_grid->extent(), _grid->_eps)) return false;

                // get the node containing the current location;
                _node = _grid->_root->leaf(r());

                // if the photon packet started outside the grid, return the corresponding nonzero-length segment;
                // otherwise fall through to determine the first actual segment
                if (ds() > 0.) return true;
            }

            // intentionally falls through
            case State::Inside:
            {
                // determine the segment from the current position to the first cell wall
                // and adjust the position and cell indices accordingly
                double xnext = (kx() < 0.0) ? _node->xmin() : _node->xmax();
                double ynext = (ky() < 0.0) ? _node->ymin() : _node->ymax();
                double znext = (kz() < 0.0) ? _node->zmin() : _node->zmax();
                double dsx = (fabs(kx()) > 1e-15) ? (xnext - rx()) / kx() : DBL_MAX;
                double dsy = (fabs(ky()) > 1e-15) ? (ynext - ry()) / ky() : DBL_MAX;
                double dsz = (fabs(kz()) > 1e-15) ? (znext - rz()) / kz() : DBL_MAX;

                double ds;
                Node::Wall wall;
                if (dsx <= dsy && dsx <= dsz)
                {
                    ds = dsx;
                    wall = (kx() < 0.0) ? Node::BACK : Node::FRONT;
                }
                else if (dsy <= dsx && dsy <= dsz)
                {
                    ds = dsy;
                    wall = (ky() < 0.0) ? Node::LEFT : Node::RIGHT;
                }
                else
                {
                    ds = dsz;
                    wall = (kz() < 0.0) ? Node::BOTTOM : Node::TOP;
                }
                propagater(ds + _grid->_eps);
                setSegment(_node->cellIndex(), ds);

                // try the most likely neighbor of the current node, and use top-down search as a fall-back
                const Node* oldnode = _node;
                _node = _node->neighbor(wall, r());
                if (!_node) _node = _grid->_root->leaf(r());

                // if we're stuck in the same node,
                // try to escape by advancing the position to the next representable coordinates
                if (_node == oldnode)
                {
                    // try to escape by advancing the position to the next representable coordinates
                    propagateToNextAfter();
                    _node = _grid->_root->leaf(r());
                }

                // if we're outside the domain or still stuck in the same node, terminate the path
                if (!_node || _node == oldnode) setState(State::Outside);
                return true;
            }

            case State::Outside:
            {
            }
        }
        return false;
    }
};

////////////////////////////////////////////////////////////////////

void AdaptiveMeshSnapshot::getEntities(EntityCollection& entities, Position bfr, Direction bfk) const
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

std::unique_ptr<PathSegmentGenerator> AdaptiveMeshSnapshot::createPathSegmentGenerator() const
{
    return std::make_unique<MySegmentGenerator>(this);
}

////////////////////////////////////////////////////////////////////
