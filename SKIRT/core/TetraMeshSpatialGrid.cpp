/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "TetraMeshSpatialGrid.hpp"
#include "Configuration.hpp"
#include "FatalError.hpp"
#include "MediumSystem.hpp"
#include "NR.hpp"
#include "Random.hpp"
#include "SiteListInterface.hpp"
#include "SpatialGridPlotFile.hpp"
#include "StringUtils.hpp"
#include "tetgen.h"

//////////////////////////////////////////////////////////////////////

namespace
{
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

    // sample random positions from the given media components with the given relative weigths
    vector<Vec> sampleMedia(const vector<Medium*>& media, const vector<double>& weights, const Box& extent,
                            int numSites)
    {
        if (media.empty()) throw FATALERROR("There is no material of the requested type in the medium system");

        // build the cumulative weight distribution
        Array Xv;
        NR::cdf(Xv, media.size(), [weights](int h) { return weights[h]; });
        auto random = media[0]->find<Random>();

        // sample
        vector<Vec> rv(numSites);
        for (int m = 0; m != numSites;)
        {
            int h = NR::locateClip(Xv, random->uniform());
            Position p = media[h]->generatePosition();
            if (extent.contains(p)) rv[m++] = p;  // discard any points outside of the domain
        }
        return rv;
    }

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

    std::array<int, 3> clockwiseVertices(int face)
    {
        std::array<int, 3> cv = {(face + 3) % 4, (face + 2) % 4, (face + 1) % 4};
        // if face is even we should swap two edges
        if (face % 2 == 0) std::swap(cv[0], cv[2]);
        return cv;
    }
}

//////////////////////////////////////////////////////////////////////

struct TetraMeshSpatialGrid::Face
{
    Face() {}

    // normals are calculated in the constructor of Tetra
    Face(int ntetra, int nface) : _ntetra(ntetra), _nface(nface) {}

    // setters & getters required?

    Vec _normal;  // outward facing normal
    int _ntetra;  // index of neighbouring tetrahedron
    int _nface;   // neighbouring face index
};

////////////////////////////////////////////////////////////////////

class TetraMeshSpatialGrid::Tetra : public Box
{
private:
    Vec _centroid;
    std::array<Vec*, 4> _vertices;
    std::array<Face, 4> _faces;

public:
    ~Tetra() {}  // vertices are deleted in ~TetraMeshSpatialGrid

    Tetra(const std::array<Vec*, 4>& vertices, const std::array<Face, 4>& faces) : _vertices(vertices), _faces(faces)
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

        // average position of all vertices
        for (int i = 0; i < 4; i++) _centroid += *_vertices[i];
        _centroid /= 4;

        // calculate normal facing out
        for (int f = 0; f < 4; f++)
        {
            std::array<int, 3> cv = clockwiseVertices(f);
            Vec e12 = edge(cv[0], cv[1]);
            Vec e13 = edge(cv[0], cv[2]);
            Vec normal = Vec::cross(e12, e13);
            normal /= normal.norm();

            Face& face = _faces[f];
            face._normal = normal;
        }
    }

    std::pair<double, int> traverse(Vec pos, Direction dir, int& enteringFace) const
    {
        int leavingFace = -1;
        double ds = DBL_MAX;

        // find entering face using a single Plücker product
        if (enteringFace == -1) enteringFace = findEnteringFace(pos, dir);

        // the translated Plücker moment in the local coordinate system
        Vec moment = Vec::cross(dir, pos - *_vertices[enteringFace]);

        // clockwise vertices around vertex 0
        std::array<int, 3> cv = clockwiseVertices(enteringFace);

        // determine orientations for use in the decision tree
        double prod0 = Vec::dot(moment, edge(cv[0], enteringFace));
        int clock0 = prod0 < 0;
        // if clockwise move clockwise else move cclockwise
        int i = clock0 ? 1 : 2;
        double prodi = Vec::dot(moment, edge(cv[i], enteringFace));
        int cclocki = prodi >= 0;

        // use plane intersection algorithm if Plücker products are ambiguous
        // this is actually more strict than the algorithm described by Maria (2017)
        // but these edge cases are incredibly rare and can cause issues
        if (prod0 == 0. || prodi == 0.)
        {
            for (int face : cv)
            {
                const Vec& n = _faces[face]._normal;
                double ndotk = Vec::dot(n, dir);
                if (ndotk > 0)
                {
                    const Vec& v = *_vertices[enteringFace];
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
            static constexpr int dtable[2][2] = {{1, 0}, {0, 2}};
            leavingFace = cv[dtable[clock0][cclocki]];
            const Vec& n = _faces[leavingFace]._normal;
            const Vec& v = *_vertices[enteringFace];
            double ndotk = Vec::dot(n, dir);
            ds = Vec::dot(n, v - pos) / ndotk;
        }

        return std::make_pair(ds, leavingFace);
    }

    ////////////////////////////////////////////////////////////////////

    int findEnteringFace(const Vec& pos, const Direction& dir) const
    {
        int enteringFace = -1;
        // clockwise and cclockwise adjacent faces when checking edge v1->v2
        static constexpr int etable[6][2] = {{3, 2}, {1, 3}, {2, 1}, {3, 0}, {0, 2}, {1, 0}};  // make this static?
        // try all 6 edges because of very rare edge cases where ray is inside edge
        // having only 1 non-zero Plücker product
        int e = 0;
        for (int v1 = 0; v1 < 3; v1++)
        {
            for (int v2 = v1 + 1; v2 < 4; v2++)
            {

                Vec moment12 = Vec::cross(dir, pos - *_vertices[v1]);
                double prod12 = Vec::dot(moment12, edge(v1, v2));
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

    ////////////////////////////////////////////////////////////////////

    bool inside(const Position& bfr) const
    {
        // since we use the k-d tree this will only slow the CellIndex
        // if (!Box::contains(bfr)) return false;

        // could optimize this slightly by using same vertex for 3 faces and do final face seperately
        // but this function acts only as a backup
        for (int f = 0; f < 4; f++)
        {
            const Face& face = _faces[f];
            const Vec* vertex = _vertices[(f + 1) % 4];  // any vertex that is on the face

            // if point->face is opposite direction as the outward pointing normal, the point is outside
            if (Vec::dot(*vertex - bfr, face._normal) < 0) return false;
        }
        return true;
    }

    double generateBarycentric(double& s, double& t, double& u) const
    {
        // https://vcg.isti.cnr.it/activities/OLD/geometryegraphics/pointintetraedro.html
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
        return 1 - u - t - s;
    }

    ////////////////////////////////////////////////////////////////////

    Position generatePosition(Random* random) const
    {
        double s = random->uniform();
        double t = random->uniform();
        double u = random->uniform();

        double r = generateBarycentric(s, t, u);

        return Position(r * *_vertices[0] + u * *_vertices[1] + t * *_vertices[2] + s * *_vertices[3]);
    }

    ////////////////////////////////////////////////////////////////////

    Vec edge(int t1, int t2) const { return *_vertices[t2] - *_vertices[t1]; }

    ////////////////////////////////////////////////////////////////////

    double volume() const { return 1 / 6. * abs(Vec::dot(Vec::cross(edge(0, 1), edge(0, 2)), edge(0, 3))); }

    double diagonal() const
    {
        double sum = 0.0;
        for (int i = 0; i < 3; ++i)
        {
            for (int j = i + 1; j < 4; ++j)
            {
                sum += edge(i, j).norm2();
            }
        }
        return sqrt(sum / 6.0);
    }

    ////////////////////////////////////////////////////////////////////

    const Face& face(int f) const { return _faces[f]; }

    ////////////////////////////////////////////////////////////////////

    const Vec& centroid() const { return _centroid; }

    ////////////////////////////////////////////////////////////////////

    const Vec& vertex(int v) const { return *_vertices[v]; }

    ////////////////////////////////////////////////////////////////////
};

//////////////////////////////////////////////////////////////////////

class TetraMeshSpatialGrid::Node
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
    Node* child(Vec bfr, const vector<const Vec*>& points) const
    {
        return lessthan(bfr, *points[_m], _axis) ? _left : _right;
    }

    // returns the other child than the one that would be apropriate for the specified query point
    Node* otherChild(Vec bfr, const vector<const Vec*>& points) const
    {
        return lessthan(bfr, *points[_m], _axis) ? _right : _left;
    }

    // returns the squared distance from the query point to the split plane
    double squaredDistanceToSplitPlane(Vec bfr, const vector<const Vec*>& points) const
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
    Node* nearest(Vec bfr, const vector<const Vec*>& points)
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

TetraMeshSpatialGrid::~TetraMeshSpatialGrid()
{
    for (auto vertex : _vertices) delete vertex;
    for (auto centroid : _centroids) delete centroid;
    for (auto tetra : _tetrahedra) delete tetra;  // vertices are already deleted
    for (auto tree : _blocktrees) delete tree;
}

//////////////////////////////////////////////////////////////////////

void TetraMeshSpatialGrid::setupSelfBefore()
{
    BoxSpatialGrid::setupSelfBefore();

    _log = find<Log>();
    _eps = 1e-12 * widths().norm();

    // determine an appropriate set of sites and construct the Tetra mesh
    switch (_policy)
    {
        case Policy::Uniform:
        {
            auto random = find<Random>();
            _vertices.resize(_numSites);
            for (int m = 0; m != _numSites; ++m) _vertices[m] = new Vec((random->position(extent())));
            break;
        }
        case Policy::CentralPeak:
        {
            auto random = find<Random>();
            const int a = 1000;  // steepness of the peak; the central 1/a portion is NOT covered
            const double rscale = extent().rmax().norm();
            _vertices.resize(_numSites);
            _vertices[0] = new Vec(0, 0, 0);
            for (int m = 1; m != _numSites;)  // skip first particle so that it remains (0,0,0)
            {
                double r = rscale * pow(1. / a, random->uniform());  // random distribution according to 1/x
                Direction k = random->direction();
                Vec p = Position(r, k);
                if (extent().contains(p)) _vertices[m++] = new Vec(p);  // discard any points outside of the domain
            }
            break;
        }
        case Policy::DustDensity:
        {
            // build a list of media that have this material type with corresponding weights
            vector<Medium*> media;
            vector<double> weights;
            auto ms = find<MediumSystem>();
            for (auto medium : ms->media())
                if (medium->mix()->isDust()) media.push_back(medium);
            for (auto medium : media) weights.push_back(medium->mass());
            vector<Vec> sites = sampleMedia(media, weights, extent(), _numSites);
            int n = sites.size();
            _vertices.resize(n);
            for (int m = 0; m != n; ++m) _vertices[m] = new Vec(sites[m]);
            break;
        }
        case Policy::ElectronDensity:
        {
            // build a list of media that have this material type with corresponding weights
            vector<Medium*> media;
            vector<double> weights;
            auto ms = find<MediumSystem>();
            for (auto medium : ms->media())
                if (medium->mix()->isElectrons()) media.push_back(medium);
            for (auto medium : media) weights.push_back(medium->number());
            vector<Vec> sites = sampleMedia(media, weights, extent(), _numSites);
            int n = sites.size();
            _vertices.resize(n);
            for (int m = 0; m != n; ++m) _vertices[m] = new Vec(sites[m]);
            break;
        }
        case Policy::GasDensity:
        {
            // build a list of media that have this material type with corresponding weights
            vector<Medium*> media;
            vector<double> weights;
            auto ms = find<MediumSystem>();
            for (auto medium : ms->media())
                if (medium->mix()->isGas()) media.push_back(medium);
            for (auto medium : media) weights.push_back(medium->number());
            vector<Vec> sites = sampleMedia(media, weights, extent(), _numSites);
            int n = sites.size();
            _vertices.resize(n);
            for (int m = 0; m != n; ++m) _vertices[m] = new Vec(sites[m]);
            break;
        }
        case Policy::ImportedSites:
        {
            auto sli = find<MediumSystem>()->interface<SiteListInterface>(2);
            // prepare the data
            int n = sli->numSites();
            _vertices.resize(n);
            for (int m = 0; m != n; ++m) _vertices[m] = new Vec(sli->sitePosition(m));
            break;
        }
    }

    if (_vertices.empty())
    {
        throw FATALERROR("No vertices available for mesh generation");
    }

    buildMesh();
    buildSearchPerBlock();
}

////////////////////////////////////////////////////////////////////

void TetraMeshSpatialGrid::buildMesh()
{
    tetgenio delaunay, refined;
    buildDelaunay(delaunay);
    if (refine())
    {
        refineDelaunay(delaunay, refined);
        storeTetrahedra(refined, true);
    }
    else
    {
        storeTetrahedra(delaunay, false);
    }
}

void TetraMeshSpatialGrid::buildDelaunay(tetgenio& out)
{
    // remove sites outside of the domain
    int numOutside = 0;
    _numVertices = _vertices.size();
    for (int m = 0; m != _numVertices; ++m)
    {
        if (!contains(*_vertices[m]))
        {
            delete _vertices[m];
            _vertices[m] = 0;
            numOutside++;
        }
    }
    if (numOutside) _numVertices = eraseNullPointers(_vertices);
    _log->info("removed " + StringUtils::toString(numOutside, 'd') + " vertices outside of the domain");

    tetgenio in;
    tetgenbehavior behavior;

    // in.firstnumber = 0; // remove me if no error
    in.numberofpoints = _vertices.size();
    in.pointlist = new REAL[in.numberofpoints * 3];
    for (int i = 0; i < in.numberofpoints; i++)
    {
        in.pointlist[i * 3 + 0] = _vertices[i]->x();
        in.pointlist[i * 3 + 1] = _vertices[i]->y();
        in.pointlist[i * 3 + 2] = _vertices[i]->z();
    }

    behavior.psc = 1;  // -s build Delaunay tetrahedralisation

    _log->info("Building Delaunay triangulation using input vertices...");
    tetrahedralize(&behavior, &in, &out);
    _log->info("Built Delaunay triangulation");
}

void TetraMeshSpatialGrid::refineDelaunay(tetgenio& in, tetgenio& out)
{
    tetgenbehavior behavior;

    // tetgen refine options
    behavior.refine = 1;   // -r
    behavior.quality = 1;  // -q
    // use default tetgen options for quality
    // correct output options for out
    behavior.neighout = 2;   // -nn
    behavior.facesout = 1;   // -f
    behavior.zeroindex = 1;  // -z

    _log->info("Refining triangulation...");
    tetrahedralize(&behavior, &in, &out);
    _log->info("Refined triangulation");
}

void TetraMeshSpatialGrid::storeTetrahedra(const tetgenio& out, bool storeVertices)
{
    // tranfser TetGen data to TetraMeshSpatialGrid data containers
    _numTetra = out.numberoftetrahedra;

    if (storeVertices)
    {
        // delete old vertices
        for (int i = 0; i < _numVertices; i++) delete _vertices[i];

        _numVertices = out.numberofpoints;

        _vertices.resize(_numVertices);
        for (int i = 0; i < _numVertices; i++)
        {
            double x = out.pointlist[3 * i + 0];
            double y = out.pointlist[3 * i + 1];
            double z = out.pointlist[3 * i + 2];

            _vertices[i] = new Vec(x, y, z);
        }
    }

    _tetrahedra.resize(_numTetra);
    for (int i = 0; i < _numTetra; i++)
    {
        std::array<Vec*, 4> vertices;
        std::array<Face, 4> faces;

        // vertices
        for (int c = 0; c < 4; c++)
        {
            vertices[c] = _vertices[out.tetrahedronlist[4 * i + c]];
        }

        // faces
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

        _centroids.push_back(&_tetrahedra[i]->centroid());
    }

    // compile statistics
    double minVol = DBL_MAX;
    double maxVol = 0.;
    double totalVol2 = 0.;
    for (int m = 0; m < _numTetra; m++)
    {
        double vol = _tetrahedra[m]->volume();
        totalVol2 += vol * vol;
        minVol = min(minVol, vol);
        maxVol = max(maxVol, vol);
    }
    double V = Box::volume();
    minVol /= V;
    maxVol /= V;
    double avgVol = 1 / (double)_numTetra;
    double varVol = (totalVol2 / _numTetra / (V * V) - avgVol * avgVol);

    // log neighbor statistics
    _log->info("Done computing tetrahedralisation");
    _log->info("  Number of vertices " + std::to_string(_numVertices));
    _log->info("  Number of tetrahedra " + std::to_string(_numTetra));
    _log->info("  Average volume fraction per cell: " + StringUtils::toString(avgVol, 'e'));
    _log->info("  Variance of volume fraction per cell: " + StringUtils::toString(varVol, 'e'));
    _log->info("  Minimum volume fraction cell: " + StringUtils::toString(minVol, 'e'));
    _log->info("  Maximum volume fraction cell: " + StringUtils::toString(maxVol, 'e'));
}

//////////////////////////////////////////////////////////////////////

TetraMeshSpatialGrid::Node* TetraMeshSpatialGrid::buildTree(vector<int>::iterator first, vector<int>::iterator last,
                                                            int depth) const
{
    auto length = last - first;
    if (length > 0)
    {
        auto median = length >> 1;
        std::nth_element(first, first + median, last, [this, depth](int m1, int m2) {
            return m1 != m2 && lessthan(*_centroids[m1], *_centroids[m2], depth % 3);
        });
        return new TetraMeshSpatialGrid::Node(*(first + median), depth, buildTree(first, first + median, depth + 1),
                                              buildTree(first + median + 1, last, depth + 1));
    }
    return nullptr;
}

////////////////////////////////////////////////////////////////////

void TetraMeshSpatialGrid::buildSearchPerBlock()
{
    // abort if there are no cells
    if (!_numTetra) return;

    _log->info("Building data structures to accelerate searching the tetrahedralisation");

    // -------------  block lists  -------------
    _nb = max(3, min(250, static_cast<int>(cbrt(_numTetra))));
    _nb2 = _nb * _nb;
    _nb3 = _nb * _nb * _nb;

    // initialize a vector of nb * nb * nb lists
    _blocklists.resize(_nb3);

    // we add the tetrahedra to all blocks they potentially overlap with
    // this will slow down the search tree but if no search tree is present
    // we can simply loop over all tetrahedra inside the block
    int i1, j1, k1, i2, j2, k2;
    for (int c = 0; c != _numTetra; ++c)
    {
        cellIndices(i1, j1, k1, _tetrahedra[c]->rmin() - Vec(_eps, _eps, _eps), _nb, _nb, _nb);
        cellIndices(i2, j2, k2, _tetrahedra[c]->rmax() + Vec(_eps, _eps, _eps), _nb, _nb, _nb);
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
    _log->info("  Number of blocks in search grid: " + std::to_string(_nb3) + " (" + std::to_string(_nb) + "^3)");
    _log->info("  Average number of cells per block: " + StringUtils::toString(avgRefsPerBlock, 'f', 1));
    _log->info("  Minimum number of cells per block: " + std::to_string(minRefsPerBlock));
    _log->info("  Maximum number of cells per block: " + std::to_string(maxRefsPerBlock));

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
    _log->info("  Number of search trees: " + std::to_string(numTrees) + " ("
               + StringUtils::toString(100. * numTrees / _nb3, 'f', 1) + "% of blocks)");
}

//////////////////////////////////////////////////////////////////////

int TetraMeshSpatialGrid::numCells() const
{
    return _numTetra;
}

//////////////////////////////////////////////////////////////////////

double TetraMeshSpatialGrid::volume(int m) const
{
    return _tetrahedra[m]->volume();
}

//////////////////////////////////////////////////////////////////////

double TetraMeshSpatialGrid::diagonal(int m) const
{
    // can use for loop but this is more readable
    const Tetra* tetra = _tetrahedra[m];
    double a = tetra->edge(0, 1).norm2();
    double b = tetra->edge(0, 2).norm2();
    double c = tetra->edge(0, 3).norm2();
    double d = tetra->edge(1, 2).norm2();
    double e = tetra->edge(1, 3).norm2();
    double f = tetra->edge(2, 3).norm2();
    return sqrt(a + b + c + d + e + f) / sqrt(6.0);
}
//////////////////////////////////////////////////////////////////////

int TetraMeshSpatialGrid::cellIndex(Position bfr) const
{
    // Ensure the position is inside the domain
    if (!contains(bfr)) return -1;

    // Determine the block in which the point falls
    int i, j, k;
    cellIndices(i, j, k, bfr, _nb, _nb, _nb);
    int b = i * _nb2 + j * _nb + k;

    // Look for the closest centroid in this block using the search tree
    Node* tree = _blocktrees[b];
    if (!tree) throw FATALERROR("No search tree found for block " + std::to_string(b));
    int m = tree->nearest(bfr, _centroids)->m();
    const Tetra* tetra = _tetrahedra[m];

    // traverse from centroid towards bfr until we pass it
    Vec pos = tetra->centroid();
    Direction dir(bfr - pos);
    double dist = dir.norm();
    dir /= dist;

    double ds;
    int leavingFace;
    int enteringFace = -1;
    while (true)
    {
        std::tie(ds, leavingFace) = tetra->traverse(pos, dir, enteringFace);

        // If no exit point was found, break
        if (leavingFace == -1) break;

        // Move to the exit point
        pos += ds * dir;
        dist -= ds;
        if (dist <= 0) return m;

        // Get neighbour information
        m = tetra->face(leavingFace)._ntetra;
        enteringFace = tetra->face(leavingFace)._nface;

        // If no next tetrahedron, break
        if (m < 0) break;

        tetra = _tetrahedra[m];
    }

    // If traversal algorithm failed to find correct cell, loop over all tetrahedra in the block
    for (int t : _blocklists[b])
        if (_tetrahedra[t]->inside(bfr)) return t;

    throw FATALERROR("Can't find random position in cell");
}

//////////////////////////////////////////////////////////////////////

Position TetraMeshSpatialGrid::centralPositionInCell(int m) const
{
    return Position(_tetrahedra[m]->centroid());
}

//////////////////////////////////////////////////////////////////////

Position TetraMeshSpatialGrid::randomPositionInCell(int m) const
{
    return _tetrahedra[m]->generatePosition(random());
}

//////////////////////////////////////////////////////////////////////

void TetraMeshSpatialGrid::writeGridPlotFiles(const SimulationItem* probe) const
{
    // create the plot files
    SpatialGridPlotFile plotxy(probe, probe->itemName() + "_grid_xy");
    SpatialGridPlotFile plotxz(probe, probe->itemName() + "_grid_xz");
    SpatialGridPlotFile plotyz(probe, probe->itemName() + "_grid_yz");
    SpatialGridPlotFile plotxyz(probe, probe->itemName() + "_grid_xyz");

    // for each site, compute the corresponding cell and output its edges
    _log->info("Writing plot files for tetrahedralisation with " + std::to_string(_numTetra) + " tetrahedra");
    _log->infoSetElapsed(_numTetra);
    int numDone = 0;
    for (int i = 0; i < _numTetra; i++)
    {
        const Tetra* tetra = _tetrahedra[i];
        vector<double> coords;
        coords.reserve(12);
        vector<int> indices;
        indices.reserve(16);

        for (int v = 0; v < 4; v++)
        {
            const Vec vertex = tetra->vertex(v);
            coords.push_back(vertex.x());
            coords.push_back(vertex.y());
            coords.push_back(vertex.z());

            // get vertices of opposite face
            std::array<int, 3> faceIndices = clockwiseVertices(v);
            indices.push_back(3);  // amount of vertices per face
            indices.push_back(faceIndices[0]);
            indices.push_back(faceIndices[1]);
            indices.push_back(faceIndices[2]);
        }

        if (tetra->zmin() <= 0 && tetra->zmax() >= 0) plotxy.writePolyhedron(coords, indices);
        if (tetra->ymin() <= 0 && tetra->ymax() >= 0) plotxz.writePolyhedron(coords, indices);
        if (tetra->xmin() <= 0 && tetra->xmax() >= 0) plotyz.writePolyhedron(coords, indices);
        if (i <= 1000)
            plotxyz.writePolyhedron(coords, indices);  // like TetraMeshSpatialGrid, but why even write at all?

        // log message if the minimum time has elapsed
        numDone++;
        if (numDone % 2000 == 0) _log->infoIfElapsed("Computed tetrehedra: ", 2000);
    }
}

//////////////////////////////////////////////////////////////////////

class TetraMeshSpatialGrid::MySegmentGenerator : public PathSegmentGenerator
{
    const TetraMeshSpatialGrid* _grid{nullptr};
    int _mr;
    int _enteringFace;

public:
    MySegmentGenerator(const TetraMeshSpatialGrid* grid) : _grid(grid) {}

    bool next() override
    {
        if (state() == State::Unknown)
        {
            // try moving the photon packet inside the grid; if this is impossible, return an empty path
            // this also changes the _state
            if (!moveInside(_grid->extent(), _grid->_eps)) return false;

            // get the index of the cell containing the current position
            _mr = _grid->cellIndex(r());
            _enteringFace = -1;

            // if the photon packet started outside the grid, return the corresponding nonzero-length segment;
            // otherwise fall through to determine the first actual segment
            if (ds() > 0.) return true;
        }

        // intentionally falls through
        if (state() == State::Inside)
        {
            double ds;
            int leavingFace;
            // loop in case no exit point was found (which should happen only rarely)
            while (true)
            {
                const Tetra* tetra = _grid->_tetrahedra[_mr];
                Position pos = r();
                Direction dir = k();

                std::tie(ds, leavingFace) = tetra->traverse(pos, dir, _enteringFace);

                // if no exit point was found, advance the current point by a small distance,
                // recalculate the cell index, and return to the start of the loop
                if (leavingFace == -1 || ds < _grid->_eps)
                {
                    propagater(_grid->_eps);
                    _mr = _grid->cellIndex(r());

                    if (_mr < 0)
                    {
                        setState(State::Outside);
                        return false;
                    }
                    else
                    {
                        _enteringFace = -1;
                    }
                }
                // otherwise set the current point to the exit point and return the path segment
                else
                {
                    propagater(ds);
                    setSegment(_mr, ds);
                    _mr = tetra->face(leavingFace)._ntetra;

                    if (_mr < 0)
                    {
                        setState(State::Outside);
                        return false;
                    }
                    else
                    {
                        _enteringFace = tetra->face(leavingFace)._nface;
                        return true;
                    }
                }
            }
        }

        if (state() == State::Outside)
        {
        }

        return false;
    }
};

//////////////////////////////////////////////////////////////////////

std::unique_ptr<PathSegmentGenerator> TetraMeshSpatialGrid::createPathSegmentGenerator() const
{
    return std::make_unique<MySegmentGenerator>(this);
}

//////////////////////////////////////////////////////////////////////
