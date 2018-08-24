/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "AdaptiveMeshSnapshot.hpp"
#include "FatalError.hpp"
#include "Log.hpp"
#include "NR.hpp"
#include "Random.hpp"
#include "SpatialGridPath.hpp"
#include "StringUtils.hpp"
#include "TextInFile.hpp"
#include "Units.hpp"

////////////////////////////////////////////////////////////////////

namespace AdaptiveMesh_Private
{
    /* The Node class is a helper class used to represent individual nodes in a tree data
       structure. An Node instance can represent a leaf or a nonleaf node. A nonleaf node
       maintains a list of pointers to its children. A leaf node instead keeps a list of pointers
       to the most likely neighbor for each of the six cell walls. */
    class Node : public Box
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
                _nodes.resize(_Nx*_Ny*_Nz);
                int m = 0;
                for (int k=0; k<_Nz; k++)
                    for (int j=0; j<_Ny; j++)
                        for (int i=0; i<_Nx; i++)
                        {
                            Vec r0 = extent.fracPos(i,j,k, _Nx, _Ny, _Nz);
                            Vec r1 = extent.fracPos(i+1,j+1,k+1, _Nx, _Ny, _Nz);
                            _nodes[m++] = new Node(Box(r0, r1), infile, leafnodes);
                        }
            }

            // if this is not a nonleaf line, it should be a leaf line
            else
            {
                _Nx = 0; _Ny = 0; _Nz = 0;

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
                _nodes[BACK]   = root->leaf(ctr + Vec(-eps, 0, 0));
                _nodes[FRONT]  = root->leaf(ctr + Vec(+eps, 0, 0));
                _nodes[LEFT]   = root->leaf(ctr + Vec(0, -eps, 0));
                _nodes[RIGHT]  = root->leaf(ctr + Vec(0, +eps, 0));
                _nodes[BOTTOM] = root->leaf(ctr + Vec(0, 0, -eps));
                _nodes[TOP]    = root->leaf(ctr + Vec(0, 0, +eps));
            }
        }

        /* The destructor releases the node's children if it is a nonleaf node. */
        ~Node()
        {
            if (!isLeaf()) for (auto node : _nodes) delete node;
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
            int i,j,k;
            cellIndices(i,j,k, r, _Nx,_Ny,_Nz);

            // get the estimated node using local Morton order
            const Node* node = _nodes[ (k*_Ny+j)*_Nx+i ];

            // if the point is NOT in the node, correct the indices and get the new node
            if (!node->contains(r))
            {
                if (r.x() < node->xmin()) i--; else if (r.x() > node->xmax()) i++;
                if (r.y() < node->ymin()) j--; else if (r.y() > node->ymax()) j++;
                if (r.z() < node->zmin()) k--; else if (r.z() > node->zmax()) k++;
                node = _nodes[ (k*_Ny+j)*_Nx+i ];
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
        enum Wall { BACK=0, FRONT, LEFT, RIGHT, BOTTOM, TOP };

        /* This function returns a pointer to the leaf node just beyond a given wall that contains
           the specified position, or null if the specified position is not inside the most likely
           neighbor for that wall, or if neighbors have not been added for the node, or if this is
           not a leaf node. */
        const Node* neighbor(Wall wall, Vec r) const
        {
            const Node* node = _nodes.size()==6 ? _nodes[wall] : nullptr;
            if (node && node->contains(r)) return node;
            else return nullptr;
        }

        /* This function returns the user properties for leaf nodes; for nonleaf nodes the
           returned array is empty. */
        const Array& properties() { return _properties; }

    private:
        int _Nx, _Ny, _Nz;   // number of grid cells in each direction; zero for leaf nodes
        int _m;              // Morton order index for the cell represented by this leaf node; -1 for nonleaf nodes
        vector<const Node*> _nodes;  // pointers to children (nonleaf nodes) or neighbors (leaf nodes)
        Array _properties;           // user-defined properties
    };
}

using namespace AdaptiveMesh_Private;

////////////////////////////////////////////////////////////////////

AdaptiveMeshSnapshot::AdaptiveMeshSnapshot()
{
}

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
        double maxT = hasTemperatureCutoff() ? maxTemperature() : 0.;

        // initialize statistics
        double totalOriginalMass = 0;
        double totalMetallicMass = 0;
        double totalEffectiveMass = 0;

        // loop over all leaf cells
        int numIgnored = 0;
        for (size_t m=0; m!=n; ++m)
        {
            const Array& prop = _cells[m]->properties();

            // original mass is zero if temperature is above cutoff or if imported mass/density is not positive
            double originalMass = 0.;
            if (maxT && prop[temperatureIndex()] > maxT) numIgnored++;
            else originalMass = max(0., massIndex()>=0 ? prop[massIndex()] : prop[densityIndex()]*_cells[m]->volume());

            double metallicMass = originalMass * (metallicityIndex()>=0 ? prop[metallicityIndex()] : 1.);
            double effectiveMass = metallicMass * multiplier();

            Mv[m] = effectiveMass;
            _rhov[m] = effectiveMass / _cells[m]->volume();

            totalOriginalMass += originalMass;
            totalMetallicMass += metallicMass;
            totalEffectiveMass += effectiveMass;
        }

        // log mass statistics
        if (numIgnored)
            log()->info("  Ignored mass in " + std::to_string(numIgnored) + " high-temperature cells" );
        log()->info("  Total original mass : " +
                    StringUtils::toString(units()->omass(totalOriginalMass),'e',4) + " " + units()->umass());
        log()->info("  Total metallic mass : "
                    + StringUtils::toString(units()->omass(totalMetallicMass),'e',4) + " " + units()->umass());
        log()->info("  Total effective mass: "
                    + StringUtils::toString(units()->omass(totalEffectiveMass),'e',4) + " " + units()->umass());

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

Box AdaptiveMeshSnapshot::extent(int m) const
{
    return _cells[m]->extent();
}

////////////////////////////////////////////////////////////////////

Vec AdaptiveMeshSnapshot::velocity(int m) const
{
    const Array& prop = _cells[m]->properties();
    return Vec(prop[velocityIndex()+0], prop[velocityIndex()+1], prop[velocityIndex()+2]);
}

////////////////////////////////////////////////////////////////////

void AdaptiveMeshSnapshot::parameters(int m, Array& params) const
{
    int n = numParameters();
    params.resize(n);
    const Array& prop = _cells[m]->properties();
    for (int i=0; i!=n; ++i) params[i] = prop[parametersIndex()+i];
}

////////////////////////////////////////////////////////////////////

Position AdaptiveMeshSnapshot::generatePosition(int m) const
{
    return random()->position(_cells[m]->extent());
}

////////////////////////////////////////////////////////////////////

double AdaptiveMeshSnapshot::density(int m) const
{
    return _rhov[m];
}

////////////////////////////////////////////////////////////////////

double AdaptiveMeshSnapshot::mass() const
{
    return _mass;
}

////////////////////////////////////////////////////////////////////

int AdaptiveMeshSnapshot::cellIndex(Position bfr) const
{
    const Node* cell = _root->leaf(bfr);
    return cell ? cell->cellIndex() : -1;
}

////////////////////////////////////////////////////////////////////

Vec AdaptiveMeshSnapshot::velocity(Position bfr) const
{
    int m = cellIndex(bfr);
    return m>=0 ? velocity(m) : Vec();
}

////////////////////////////////////////////////////////////////////

double AdaptiveMeshSnapshot::density(Position bfr) const
{
    int m = cellIndex(bfr);
    return m>=0 ? _rhov[m] : 0;
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

void AdaptiveMeshSnapshot::path(SpatialGridPath* path) const
{
    // initialize the path
    path->clear();

    // if the photon package starts outside the dust grid, move it into the first grid cell that it will pass
    Position r = path->moveInside(_root->extent(), _eps);

    // get the node containing the current location;
    // if the position is not inside the grid, return an empty path
    const Node* node = _root->leaf(r);
    if (!node) return path->clear();

    // start the loop over nodes/path segments until we leave the grid
    double kx,ky,kz;
    path->direction().cartesian(kx,ky,kz);
    while (node)
    {
        double xnext = (kx<0.0) ? node->xmin() : node->xmax();
        double ynext = (ky<0.0) ? node->ymin() : node->ymax();
        double znext = (kz<0.0) ? node->zmin() : node->zmax();
        double dsx = (fabs(kx)>1e-15) ? (xnext-r.x())/kx : DBL_MAX;
        double dsy = (fabs(ky)>1e-15) ? (ynext-r.y())/ky : DBL_MAX;
        double dsz = (fabs(kz)>1e-15) ? (znext-r.z())/kz : DBL_MAX;

        double ds;
        Node::Wall wall;
        if (dsx<=dsy && dsx<=dsz)
        {
            ds = dsx;
            wall = (kx<0.0) ? Node::BACK : Node::FRONT;
        }
        else if (dsy<=dsx && dsy<=dsz)
        {
            ds = dsy;
            wall = (ky<0.0) ? Node::LEFT : Node::RIGHT;
        }
        else
        {
            ds = dsz;
            wall = (kz<0.0) ? Node::BOTTOM : Node::TOP;
        }
        path->addSegment(node->cellIndex(), ds);
        r += (ds+_eps)*(path->direction());

        // try the most likely neighbor of the current node, and use top-down search as a fall-back
        const Node* oldnode = node;
        node = node->neighbor(wall,r);
        if (!node) node = _root->leaf(r);

        // if we're stuck in the same node...
        if (node==oldnode)
        {
            // try to escape by advancing the position to the next representable coordinates
            r.set( nextafter(r.x(), (kx<0.0) ? -DBL_MAX : DBL_MAX),
                   nextafter(r.y(), (ky<0.0) ? -DBL_MAX : DBL_MAX),
                   nextafter(r.z(), (kz<0.0) ? -DBL_MAX : DBL_MAX) );
            node = _root->leaf(r);

            // if that didn't work, terminate the path
            if (node==oldnode)
            {
                log()->warning("Photon packet is stuck in cell "
                                + std::to_string(node->cellIndex()) + " -- terminating this path");
                break;
            }
        }
    }
}

////////////////////////////////////////////////////////////////////
