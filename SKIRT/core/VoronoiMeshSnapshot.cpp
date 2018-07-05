/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "VoronoiMeshSnapshot.hpp"
#include "Log.hpp"
#include "NR.hpp"
#include "Random.hpp"
#include "SiteListInterface.hpp"
#include "StringUtils.hpp"
#include "TextInFile.hpp"
#include "Units.hpp"
#include "container.hh"

////////////////////////////////////////////////////////////////////

namespace VoronoiMesh_Private
{
    // class to hold the information about a Voronoi cell that is relevant for calculating paths and densities
    class VoronoiCell : public Box  // enclosing box
    {
    private:
        Vec _r;                     // site position
        Vec _c;                     // centroid position
        double _volume;             // volume
        vector<int> _neighbors;     // list of neighbor indices in cells vector

    public:
        // constructor stores the specified site position; the other data members are set to zero
        VoronoiCell(Vec r) : _r(r), _volume(0) { }

        // initializes the receiver with information taken from the specified fully computed Voronoi cell
        void init(voro::voronoicell_neighbor& cell)
        {
            // copy basic geometric info
            double cx, cy, cz;
            cell.centroid(cx,cy,cz);
            _c = Vec(cx,cy,cz) + _r;
            _volume = cell.volume();

            // get the minimal and maximal coordinates of the box enclosing the cell
            vector<double> coords;
            cell.vertices(_r.x(),_r.y(),_r.z(), coords);
            double xmin = DBL_MAX;  double ymin = DBL_MAX;  double zmin = DBL_MAX;
            double xmax = -DBL_MAX; double ymax = -DBL_MAX; double zmax = -DBL_MAX;
            int n = coords.size();
            for (int i=0; i<n; i+=3)
            {
                xmin = min(xmin,coords[i]); ymin = min(ymin,coords[i+1]); zmin = min(zmin,coords[i+2]);
                xmax = max(xmax,coords[i]); ymax = max(ymax,coords[i+1]); zmax = max(zmax,coords[i+2]);
            }

            // set our inherited Box to this bounding box
            setExtent(xmin, ymin, zmin, xmax, ymax, zmax);

            // copy a list of neighboring cell/site ids
            cell.neighbors(_neighbors);
        }

        // returns the cell's site position
        Vec position() const { return _r; }

        // returns the squared distance from the cell's site to the specified point
        double squaredDistanceTo(Vec r) const { return (r-_r).norm2(); }

        // returns the central position in the cell
        Vec centroid() const { return _c; }

        // returns the volume of the cell; overriding volume() function of Box bas class
        double volume() const { return _volume; }

        // returns a list of neighboring cell/site ids
        const vector<int>& neighbors() { return _neighbors; }
    };

    // function to compare two points according to the specified axis (0,1,2)
    bool lessthan(Vec p1, Vec p2, int axis)
    {
        switch (axis)
        {
        case 0:  // split on x
            if (p1.x()<p2.x()) return true;
            if (p1.x()>p2.x()) return false;
            if (p1.y()<p2.y()) return true;
            if (p1.y()>p2.y()) return false;
            if (p1.z()<p2.z()) return true;
            return false;
        case 1:  // split on y
            if (p1.y()<p2.y()) return true;
            if (p1.y()>p2.y()) return false;
            if (p1.z()<p2.z()) return true;
            if (p1.z()>p2.z()) return false;
            if (p1.x()<p2.x()) return true;
            return false;
        case 2:  // split on z
            if (p1.z()<p2.z()) return true;
            if (p1.z()>p2.z()) return false;
            if (p1.x()<p2.x()) return true;
            if (p1.x()>p2.x()) return false;
            if (p1.y()<p2.y()) return true;
            return false;
        default: // this should never happen
            return false;
        }
    }
}

////////////////////////////////////////////////////////////////////

namespace VoronoiMesh_Private
{
    // class to hold a node in the binary search tree (see en.wikipedia.org/wiki/Kd-tree)
    class Node
    {
    private:
        int _m;         // index in _cells to the site defining the split at this node
        int _axis;      // split axis for this node (0,1,2)
        Node* _up;      // ptr to the parent node
        Node* _left;    // ptr to the left child node
        Node* _right;   // ptr to the right child node

        // returns the square of its argument
        static double sqr(double x) { return x*x; }

    public:
        // constructor stores the specified site index and child pointers (which may be null)
        Node(int m, int depth, Node* left, Node* right) : _m(m), _axis(depth%3), _up(0), _left(left), _right(right)
        {
            if (_left) _left->setParent(this);
            if (_right) _right->setParent(this);
        }

        // destructor destroys the children
        ~Node() { delete _left; delete _right; }

        // sets parent pointer (called from parent's constructor)
        void setParent(Node* up) { _up = up; }

        // returns the corresponding data member
        int m() const { return _m; }
        Node* up() const { return _up; }
        Node* left() const { return _left; }
        Node* right() const { return _right; }

        // returns the apropriate child for the specified query point
        Node* child(Vec bfr, const vector<VoronoiCell*>& cells) const
            { return lessthan(bfr, cells[_m]->position(), _axis) ? _left : _right; }

        // returns the other child than the one that would be apropriate for the specified query point
        Node* otherChild(Vec bfr, const vector<VoronoiCell*>& cells) const
            { return lessthan(bfr, cells[_m]->position(), _axis) ? _right : _left; }

        // returns the squared distance from the query point to the split plane
        double squaredDistanceToSplitPlane(Vec bfr, const vector<VoronoiCell*>& cells) const
        {
            switch (_axis)
            {
            case 0:  // split on x
                return sqr(cells[_m]->position().x() - bfr.x());
            case 1:  // split on y
                return sqr(cells[_m]->position().y() - bfr.y());
            case 2:  // split on z
                return sqr(cells[_m]->position().z() - bfr.z());
            default: // this should never happen
                return 0;
            }
        }

        // returns the node in this subtree that represents the site nearest to the query point
        Node* nearest(Vec bfr, const vector<VoronoiCell*>& cells)
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
                if (current==this) break;
                current = current->up();
            }
            return best;
        }
    };
}

using namespace VoronoiMesh_Private;

////////////////////////////////////////////////////////////////////

VoronoiMeshSnapshot::VoronoiMeshSnapshot()
{
}

////////////////////////////////////////////////////////////////////

VoronoiMeshSnapshot::~VoronoiMeshSnapshot()
{
    for (auto cell : _cells) delete cell;
    for (auto tree : _blocktrees) delete tree;
}

////////////////////////////////////////////////////////////////////

void VoronoiMeshSnapshot::open(const SimulationItem* item, string filename, string description)
{
    Snapshot::open(item, filename, description);
    importPosition();
}

////////////////////////////////////////////////////////////////////

void VoronoiMeshSnapshot::readAndClose()
{
    // read the site info into memory
    //  -> if the user configured a temperature cutoff, we need to skip the "hot" sites
    //  -> in any case, we need to skip sites outside our domain
    double maxT  = hasTemperatureCutoff() ? maxTemperature() : 0.;
    int numIgnored = 0;
    Array prop;
    while (infile()->readRow(prop))
    {
        Vec r(prop[positionIndex()+0], prop[positionIndex()+1], prop[positionIndex()+2]);
        if (!_extent.contains(r) || (maxT && prop[temperatureIndex()] > maxT)) numIgnored++;
        else _propv.push_back(prop);
    }

    // close the file
    Snapshot::readAndClose();

    // log the number of sites
    if (!numIgnored)
    {
        log()->info("  Number of sites: " + std::to_string(_propv.size()));
    }
    else
    {
        log()->info("  Number of sites ignored: " + std::to_string(numIgnored));
        log()->info("  Number of sites retained: " + std::to_string(_propv.size()));
    }

    // calculate the Voronoi cells
    vector<Vec> sites;
    sites.reserve(_propv.size());
    for (const Array& prop : _propv)
        sites.emplace_back(prop[positionIndex()+0], prop[positionIndex()+1], prop[positionIndex()+2]);
    buildMesh(sites);

    // if a mass density policy has been set, calculate masses and densities for all cells
    if (hasMassDensityPolicy())
    {
        // allocate vectors for mass and density
        size_t n = _propv.size();
        Array Mv(n);
        _rhov.resize(n);

        // initialize statistics
        double totalOriginalMass = 0;
        double totalMetallicMass = 0;
        double totalEffectiveMass = 0;

        // loop over all sites/cells
        for (size_t m=0; m!=n; ++m)
        {
            const Array& prop = _propv[m];
            double originalMass = massIndex() ? prop[massIndex()] : prop[densityIndex()] * _cells[m]->volume();
            double metallicMass = originalMass * (metallicityIndex()>=0 ? prop[metallicityIndex()] : 1.);
            double effectiveMass = metallicMass * multiplier();

            Mv[m] = effectiveMass;
            _rhov[m] = effectiveMass / _cells[m]->volume();

            totalOriginalMass += originalMass;
            totalMetallicMass += metallicMass;
            totalEffectiveMass += effectiveMass;
        }

        // log mass statistics
        log()->info("  Total original mass : " +
                    StringUtils::toString(units()->omass(totalOriginalMass),'e',4) + " " + units()->umass());
        log()->info("  Total metallic mass : "
                    + StringUtils::toString(units()->omass(totalMetallicMass),'e',4) + " " + units()->umass());
        log()->info("  Total effective mass: "
                    + StringUtils::toString(units()->omass(totalEffectiveMass),'e',4) + " " + units()->umass());

        // if one of the total masses is negative, suppress the complete mass distribution
        if (totalOriginalMass < 0 || totalMetallicMass < 0 || totalEffectiveMass < 0)
        {
            log()->warning("  Total imported mass is negative; suppressing the complete mass distribution");
            totalEffectiveMass = 0;
            _propv.clear();
            _rhov.resize(0);
            _cells.clear();
        }

        // remember the effective mass
        _mass = totalEffectiveMass;

        // if there is at least one site...
        if (!_propv.empty())
        {
            // construct a vector with the normalized cumulative site densities
            NR::cdf(_cumrhov, Mv);

            // build the search data structure
            buildSearch();
        }
    }
}

////////////////////////////////////////////////////////////////////

void VoronoiMeshSnapshot::setExtent(const Box& extent)
{
    _extent = extent;
    _eps = 1e-12 * extent.widths().norm();
}

////////////////////////////////////////////////////////////////////

VoronoiMeshSnapshot::VoronoiMeshSnapshot(const SimulationItem* item, const Box& extent, std::string filename)
{
    setContext(item);
    setExtent(extent);

    // read the input file
    TextInFile in(item, filename, "Voronoi sites");
    in.addColumn("position x", "length", "pc");
    in.addColumn("position y", "length", "pc");
    in.addColumn("position z", "length", "pc");
    auto propv = in.readAllRows();
    in.close();

    // prepare the data
    vector<Vec> sites;
    sites.reserve(_propv.size());
    for (const Array& prop : propv) sites.emplace_back(prop[0], prop[1], prop[2]);

    // calculate the Voronoi cells
    buildMesh(sites);
    buildSearch();
}

////////////////////////////////////////////////////////////////////

VoronoiMeshSnapshot::VoronoiMeshSnapshot(const SimulationItem* item, const Box& extent, SiteListInterface* sli)
{
    setContext(item);
    setExtent(extent);

    // prepare the data
    int n = sli->numSites();
    vector<Vec> sites(n);
    for (int m=0; m!=n; ++m) sites[m] = sli->sitePosition(m);

    // calculate the Voronoi cells
    buildMesh(sites);
    buildSearch();
}

////////////////////////////////////////////////////////////////////

VoronoiMeshSnapshot::VoronoiMeshSnapshot(const SimulationItem* item, const Box& extent, const vector<Vec>& sites)
{
    setContext(item);
    setExtent(extent);

    // calculate the Voronoi cells
    buildMesh(sites);
    buildSearch();
}

////////////////////////////////////////////////////////////////////

void VoronoiMeshSnapshot::buildMesh(const vector<Vec>& sites)
{

}

////////////////////////////////////////////////////////////////////

void VoronoiMeshSnapshot::buildSearch()
{

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

}

////////////////////////////////////////////////////////////////////

double VoronoiMeshSnapshot::volume(int m) const
{

}

////////////////////////////////////////////////////////////////////

Box VoronoiMeshSnapshot::extent(int m) const
{

}

////////////////////////////////////////////////////////////////////

Vec VoronoiMeshSnapshot::velocity(int m) const
{
    return Vec(_propv[m][velocityIndex()+0], _propv[m][velocityIndex()+1], _propv[m][velocityIndex()+2]);
}

////////////////////////////////////////////////////////////////////

void VoronoiMeshSnapshot::parameters(int m, Array& params) const
{
    int n = numParameters();
    params.resize(n);
    for (int i=0; i!=n; ++i) params[i] = _propv[m][parametersIndex()+i];
}

////////////////////////////////////////////////////////////////////

Position VoronoiMeshSnapshot::generatePosition(int m) const
{
/*

    // get center position and size for this particle
    Position rc(_propv[m][positionIndex()+0], _propv[m][positionIndex()+1], _propv[m][positionIndex()+2]);
    double h = _propv[m][sizeIndex()];

    // sample random position inside the smoothed unit volume
    double u = _kernel->generateRadius();
    Direction k = random()->direction();

    return Position(rc + k*u*h);
    */
}

////////////////////////////////////////////////////////////////////

double VoronoiMeshSnapshot::density(int m) const
{

}

////////////////////////////////////////////////////////////////////

double VoronoiMeshSnapshot::mass() const
{
    return _mass;
}

////////////////////////////////////////////////////////////////////

int VoronoiMeshSnapshot::cellIndex(Position bfr) const
{

}

////////////////////////////////////////////////////////////////////

double VoronoiMeshSnapshot::density(Position bfr) const
{
/*    double sum = 0.;
    if (_grid) for (const SmoothedParticle* p : _grid->particlesFor(bfr))
    {
        double u = (bfr - p->center()).norm() / p->radius();
        sum += _kernel->density(u) * p->mass();
    }
    return sum > 0. ? sum : 0.;     // guard against negative densities
*/
}

////////////////////////////////////////////////////////////////////

Position VoronoiMeshSnapshot::generatePosition() const
{
    // if there are no particles, return the origin
    if (_propv.empty()) return Position();

    // select a particle according to its mass contribution
    int m = NR::locateClip(_cumrhov, random()->uniform());

    return generatePosition(m);
}

////////////////////////////////////////////////////////////////////
