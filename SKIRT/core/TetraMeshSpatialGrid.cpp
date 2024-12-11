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

using std::array;

//////////////////////////////////////////////////////////////////////

namespace
{
    // returns the linear index for element (i,j,k) in a p*p*p table
    inline int index(int p, int i, int j, int k)
    {
        return ((i * p) + j) * p + k;
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

    inline std::array<int, 3> clockwiseVertices(int face)
    {
        std::array<int, 3> cv = {(face + 3) % 4, (face + 2) % 4, (face + 1) % 4};
        // if face is even we should swap two edges
        if (face % 2 == 0) std::swap(cv[0], cv[2]);
        return cv;
    }
}

//////////////////////////////////////////////////////////////////////

// This is a helper class for organizing cuboidal cells in a smart grid, so that
// it is easy to retrieve the first cell that overlaps a given point in space.
// The Box object on which this class is based specifies a cuboid guaranteed to
// enclose all cells in the grid.
class TetraMeshSpatialGrid::CellGrid
{
    // data members initialized during construction
    const vector<Tetra>& _tetra;   // reference to the original list of tetrahedra
    int _gridsize;                 // number of grid cells in each spatial direction
    Array _xgrid, _ygrid, _zgrid;  // the m+1 grid separation points for each spatial direction
    vector<vector<int>> _listv;    // the m*m*m lists of indices for cells overlapping each grid cell
    int _pmin, _pmax, _ptotal;     // minimum, maximum nr of cells in list; total nr of cells in listv

public:
    // The constructor creates a cuboidal grid of the specified number of grid cells in each
    // spatial direction, and for each of the grid cells it builds a list of all cells
    // overlapping the grid cell. In an attempt to distribute the cells evenly over the
    // grid cells, the sizes of the grid cells in each spatial direction are chosen so that
    // the cell centers are evenly distributed over the grid cells.
    CellGrid(const vector<Tetra>& tetra, Box extent, int gridsize) : _tetra(tetra), _gridsize(gridsize)
    {
        // build the grids in each spatial direction
        makegrid(0, gridsize, _xgrid, extent.xmin(), extent.xmax());
        makegrid(1, gridsize, _ygrid, extent.ymin(), extent.ymax());
        makegrid(2, gridsize, _zgrid, extent.zmin(), extent.zmax());

        // make room for p*p*p grid cells
        _listv.resize(gridsize * gridsize * gridsize);

        // add each cell to the list for every grid cell that it overlaps
        int n = _tetra.size();
        for (int m = 0; m != n; ++m)
        {
            Box boundingBox = _tetra[m].extent();

            // find indices for first and last grid cell overlapped by cell, in each spatial direction
            int i1 = NR::locateClip(_xgrid, boundingBox.xmin());
            int i2 = NR::locateClip(_xgrid, boundingBox.xmax());
            int j1 = NR::locateClip(_ygrid, boundingBox.ymin());
            int j2 = NR::locateClip(_ygrid, boundingBox.ymax());
            int k1 = NR::locateClip(_zgrid, boundingBox.zmin());
            int k2 = NR::locateClip(_zgrid, boundingBox.zmax());

            // add the cell to all grid cells in that 3D range
            for (int i = i1; i <= i2; i++)
                for (int j = j1; j <= j2; j++)
                    for (int k = k1; k <= k2; k++)
                    {
                        _listv[index(gridsize, i, j, k)].push_back(m);
                    }
        }

        // calculate statistics
        _pmin = n;
        _pmax = 0;
        _ptotal = 0;
        for (int index = 0; index < gridsize * gridsize * gridsize; index++)
        {
            int size = _listv[index].size();
            _pmin = min(_pmin, size);
            _pmax = max(_pmax, size);
            _ptotal += size;
        }
    }

    void makegrid(int axis, int gridsize, Array& grid, double cmin, double cmax)
    {
        int n = _tetra.size();

        // determine the cell distribution by binning at a decent resolution
        int nbins = gridsize * 100;
        double binwidth = (cmax - cmin) / nbins;
        vector<int> bins(nbins);
        for (const Tetra& tetra : _tetra)
        {
            double center = tetra.centroid(axis);
            bins[static_cast<int>((center - cmin) / binwidth)] += 1;
        }

        // determine grid separation points based on the cumulative distribution
        grid.resize(gridsize + 1);
        grid[0] = -std::numeric_limits<double>::infinity();
        int percell = n / gridsize;  // target number of particles per cell
        int cumul = 0;               // cumulative number of particles in processed bins
        int gridindex = 1;           // index of the next grid separation point to be filled
        for (int binindex = 0; binindex < nbins; binindex++)
        {
            cumul += bins[binindex];
            if (cumul > percell * gridindex)
            {
                grid[gridindex] = cmin + (binindex + 1) * binwidth;
                gridindex += 1;
                if (gridindex >= gridsize) break;
            }
        }
        grid[gridsize] = std::numeric_limits<double>::infinity();
    }

    // This function returns the smallest number of cells overlapping a single grid cell.
    int minCellRefsPerCell() const { return _pmin; }

    // This function returns the largest number of cells overlapping a single grid cell.
    int maxCellRefsPerCell() const { return _pmax; }

    // This function returns the total number of cell references for all cells in the grid.
    int totalCellRefs() const { return _ptotal; }

    // This function returns the index (in the list originally passed to the constructor)
    // of the first cell in the list that overlaps the specified position,
    // or -1 if none of the cells in the list overlap the specified position.
    int cellIndexFor(Position r) const
    {
        // locate the grid cell containing the specified position
        int i = NR::locateClip(_xgrid, r.x());
        int j = NR::locateClip(_ygrid, r.y());
        int k = NR::locateClip(_zgrid, r.z());

        // search the list of cells for that grid cell
        for (int m : _listv[index(_gridsize, i, j, k)])
        {
            if (_tetra[m].contains(r)) return m;
        }
        return -1;
    }
};

//////////////////////////////////////////////////////////////////////

struct TetraMeshSpatialGrid::Face
{
    Face() {};

    Face(int ntetra, int nface, Vec normal) : _ntetra(ntetra), _nface(nface), _normal(normal) {}

    Vec _normal;  // outward facing normal
    int _ntetra;  // index of neighbouring tetrahedron
    int _nface;   // neighbouring face index
};

////////////////////////////////////////////////////////////////////

class TetraMeshSpatialGrid::Tetra
{
private:
    const vector<Vec>& _vertices;  // reference to the full list of vertices
    Box _extent;                   // bounding box of the tetrahedron
    Vec _centroid;                 // barycenter of the tetrahedron
    array<int, 4> _vertexIndices;  // indices of the vertices in the full list
    array<Face, 4> _faces;         // face information

public:
    Tetra(const vector<Vec>& vertices, const array<int, 4>& vertexIndices, const array<Face, 4>& faces)
        : _vertices(vertices), _vertexIndices(vertexIndices), _faces(faces)
    {
        double xmin = DBL_MAX;
        double ymin = DBL_MAX;
        double zmin = DBL_MAX;
        double xmax = -DBL_MAX;
        double ymax = -DBL_MAX;
        double zmax = -DBL_MAX;
        for (int vi : _vertexIndices)
        {
            const Vec& vertex = vertices[vi];

            xmin = min(xmin, vertex.x());
            ymin = min(ymin, vertex.y());
            zmin = min(zmin, vertex.z());
            xmax = max(xmax, vertex.x());
            ymax = max(ymax, vertex.y());
            zmax = max(zmax, vertex.z());

            _centroid += vertex;
        }
        // set bounding box
        _extent = Box(xmin, ymin, zmin, xmax, ymax, zmax);

        // average position of all vertices
        _centroid /= 4;
    }

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
                Vec moment12 = Vec::cross(dir, pos - vertex(v1));
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

    bool contains(const Position& bfr) const
    {
        if (!_extent.contains(bfr)) return false;

        // could optimize this slightly by using same vertex for 3 faces and do final face seperately, but this is more readable
        for (int f = 0; f < 4; f++)
        {
            const Face& face = _faces[f];
            Vec v = vertex((f + 1) % 4);  // any vertex that is on the face

            // if point->face is opposite direction as the outward pointing normal, the point is outside
            if (Vec::dot(v - bfr, face._normal) < 0) return false;
        }
        return true;
    }

    // https://vcg.isti.cnr.it/activities/OLD/geometryegraphics/pointintetraedro.html
    double generateBarycentric(double& s, double& t, double& u) const
    {
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

    Position generatePosition(Random* random) const
    {
        double s = random->uniform();
        double t = random->uniform();
        double u = random->uniform();

        double r = generateBarycentric(s, t, u);

        return Position(r * vertex(0) + u * vertex(1) + t * vertex(2) + s * vertex(3));
    }

    // get the vertex using local indices [0, 3]
    Vec vertex(int t) const { return _vertices[_vertexIndices[t]]; }

    // get the edge t1->t2 using local indices [0, 3]
    Vec edge(int t1, int t2) const { return vertex(t2) - vertex(t1); }

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

    const array<Face, 4>& faces() const { return _faces; }

    const Vec& centroid() const { return _centroid; }

    const double centroid(int axis) const
    {
        switch (axis % 3)
        {
            case 0: return _centroid.x();
            case 1: return _centroid.y();
            case 2: return _centroid.z();
            default: return 0.0;
        }
    }

    const Box& extent() const { return _extent; }
};

////////////////////////////////////////////////////////////////////

TetraMeshSpatialGrid::~TetraMeshSpatialGrid()
{
    delete _grid;
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
            for (int m = 0; m != _numSites; ++m) _vertices[m] = random->position(extent());
            break;
        }
        case Policy::CentralPeak:
        {
            auto random = find<Random>();
            const int a = 1000;  // steepness of the peak; the central 1/a portion is NOT covered
            const double rscale = extent().rmax().norm();
            _vertices.resize(_numSites);
            _vertices[0] = Vec(0, 0, 0);
            for (int m = 1; m != _numSites;)  // skip first particle so that it remains (0,0,0)
            {
                double r = rscale * pow(1. / a, random->uniform());  // random distribution according to 1/x
                Direction k = random->direction();
                Position p(r, k);
                if (extent().contains(p)) _vertices[m++] = p;  // discard any points outside of the domain
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
            for (int m = 0; m != n; ++m) _vertices[m] = sites[m];
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
            for (int m = 0; m != n; ++m) _vertices[m] = sites[m];
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
            for (int m = 0; m != n; ++m) _vertices[m] = sites[m];
            break;
        }
        case Policy::ImportedSites:
        {
            auto sli = find<MediumSystem>()->interface<SiteListInterface>(2);
            // prepare the data
            int n = sli->numSites();
            _vertices.resize(n);
            for (int m = 0; m != n; ++m) _vertices[m] = sli->sitePosition(m);
            break;
        }
    }

    if (_vertices.empty())
    {
        throw FATALERROR("No vertices available for mesh generation");
    }

    _numVertices = _vertices.size();

    buildMesh();
    buildSearch();
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

////////////////////////////////////////////////////////////////////

void TetraMeshSpatialGrid::buildDelaunay(tetgenio& out)
{
    // remove vertices outside of the domain
    auto sitesEnd = std::remove_if(_vertices.begin(), _vertices.end(), [this](const Vec& vertex) {
        if (!contains(vertex)) return true;  // remove vertex
        return false;
    });
    _vertices.erase(sitesEnd, _vertices.end());
    int numOutside = _numVertices - _vertices.size();
    _numVertices = _vertices.size();
    if (numOutside) _log->info("removed " + StringUtils::toString(numOutside, 'd') + " vertices outside of the domain");

    tetgenio in;
    tetgenbehavior behavior;

    // in.firstnumber = 0; // remove me if no error
    in.numberofpoints = _vertices.size();
    in.pointlist = new REAL[in.numberofpoints * 3];
    for (int i = 0; i < in.numberofpoints; i++)
    {
        in.pointlist[i * 3 + 0] = _vertices[i].x();
        in.pointlist[i * 3 + 1] = _vertices[i].y();
        in.pointlist[i * 3 + 2] = _vertices[i].z();
    }

    behavior.psc = 1;  // -s build Delaunay tetrahedralisation

    _log->info("Building Delaunay triangulation using input vertices...");
    tetrahedralize(&behavior, &in, &out);
    _log->info("Built Delaunay triangulation");
}

////////////////////////////////////////////////////////////////////

void TetraMeshSpatialGrid::refineDelaunay(tetgenio& in, tetgenio& out)
{
    tetgenbehavior behavior;

    // tetgen refine options
    behavior.refine = 1;   // -r
    behavior.quality = 1;  // -q with default tetgen options for quality
    // correct output options for out
    behavior.neighout = 2;   // -nn
    behavior.facesout = 1;   // -f
    behavior.zeroindex = 1;  // -z

    _log->info("Refining triangulation...");
    tetrahedralize(&behavior, &in, &out);
    _log->info("Refined triangulation");
}

////////////////////////////////////////////////////////////////////

void TetraMeshSpatialGrid::storeTetrahedra(const tetgenio& out, bool storeVertices)
{
    // tranfser TetGen data to TetraMeshSpatialGrid data containers
    _numCells = out.numberoftetrahedra;

    // replace old vertices
    if (storeVertices)
    {
        _numVertices = out.numberofpoints;

        _vertices.resize(_numVertices);
        for (int i = 0; i < _numVertices; i++)
        {
            double x = out.pointlist[3 * i + 0];
            double y = out.pointlist[3 * i + 1];
            double z = out.pointlist[3 * i + 2];

            _vertices[i] = Vec(x, y, z);
        }
    }

    _tetrahedra.reserve(_numCells);  // can't make empty Tetra so don't use resize
    for (int i = 0; i < _numCells; i++)
    {
        std::array<int, 4> vertexIndices;
        std::array<Face, 4> faces;

        // vertices
        for (int c = 0; c < 4; c++)
        {
            vertexIndices[c] = out.tetrahedronlist[4 * i + c];
        }

        // faces
        for (int f = 0; f < 4; f++)
        {
            // -1 if no neighbor
            int ntetra = out.neighborlist[4 * i + f];

            // find which face is shared with neighbor
            int nface = -1;
            if (ntetra != -1)
            {
                for (int fn = 0; fn < 4; fn++)
                {
                    if (out.neighborlist[4 * ntetra + fn] == i)
                    {
                        nface = fn;
                        break;
                    }
                }
            }

            // compute outward facing normal of face
            std::array<int, 3> cv = clockwiseVertices(f);
            Vec v0 = _vertices[vertexIndices[cv[0]]];
            Vec e12 = _vertices[vertexIndices[cv[1]]] - v0;
            Vec e13 = _vertices[vertexIndices[cv[2]]] - v0;
            Vec normal = Vec::cross(e12, e13);
            normal /= normal.norm();

            faces[f] = Face(ntetra, nface, normal);
        }

        _tetrahedra.emplace_back(_vertices, vertexIndices, faces);
    }

    // compile statistics
    double minVol = DBL_MAX;
    double maxVol = 0.;
    double totalVol2 = 0.;
    for (int m = 0; m < _numCells; m++)
    {
        double vol = _tetrahedra[m].volume();
        totalVol2 += vol * vol;
        minVol = min(minVol, vol);
        maxVol = max(maxVol, vol);
    }
    double V = Box::volume();
    minVol /= V;
    maxVol /= V;
    double avgVol = 1 / (double)_numCells;
    double varVol = (totalVol2 / _numCells / (V * V) - avgVol * avgVol);

    // log neighbor statistics
    _log->info("Done computing tetrahedralisation");
    _log->info("  Number of vertices " + std::to_string(_numVertices));
    _log->info("  Number of tetrahedra " + std::to_string(_numCells));
    _log->info("  Average volume fraction per cell: " + StringUtils::toString(avgVol, 'e'));
    _log->info("  Variance of volume fraction per cell: " + StringUtils::toString(varVol, 'e'));
    _log->info("  Minimum volume fraction cell: " + StringUtils::toString(minVol, 'e'));
    _log->info("  Maximum volume fraction cell: " + StringUtils::toString(maxVol, 'e'));
}

////////////////////////////////////////////////////////////////////

void TetraMeshSpatialGrid::buildSearch()
{
    int gridsize = max(20, static_cast<int>(pow(_tetrahedra.size(), 1. / 3.) / 5));
    string size = std::to_string(gridsize);
    _log->info("Constructing intermediate " + size + "x" + size + "x" + size + " grid for cells...");
    _grid = new CellGrid(_tetrahedra, extent(), gridsize);
    _log->info("  Smallest number of cells per grid cell: " + std::to_string(_grid->minCellRefsPerCell()));
    _log->info("  Largest  number of cells per grid cell: " + std::to_string(_grid->maxCellRefsPerCell()));
    _log->info("  Average  number of cells per grid cell: "
               + StringUtils::toString(_grid->totalCellRefs() / double(gridsize * gridsize * gridsize), 'f', 1));
}

//////////////////////////////////////////////////////////////////////

int TetraMeshSpatialGrid::numCells() const
{
    return _numCells;
}

//////////////////////////////////////////////////////////////////////

double TetraMeshSpatialGrid::volume(int m) const
{
    return _tetrahedra[m].volume();
}

//////////////////////////////////////////////////////////////////////

double TetraMeshSpatialGrid::diagonal(int m) const
{
    return _tetrahedra[m].diagonal();
}
//////////////////////////////////////////////////////////////////////

int TetraMeshSpatialGrid::cellIndex(Position bfr) const
{
    return _grid->cellIndexFor(bfr);
}

//////////////////////////////////////////////////////////////////////

Position TetraMeshSpatialGrid::centralPositionInCell(int m) const
{
    return Position(_tetrahedra[m].centroid());
}

//////////////////////////////////////////////////////////////////////

Position TetraMeshSpatialGrid::randomPositionInCell(int m) const
{
    return _tetrahedra[m].generatePosition(random());
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
    _log->info("Writing plot files for tetrahedralisation with " + std::to_string(_numCells) + " tetrahedra");
    _log->infoSetElapsed(_numCells);
    int numDone = 0;
    for (int i = 0; i < _numCells; i++)
    {
        const Tetra& tetra = _tetrahedra[i];
        vector<double> coords;
        coords.reserve(12);
        vector<int> indices;
        indices.reserve(16);

        for (int v = 0; v < 4; v++)
        {
            const Vec vertex = tetra.vertex(v);
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

        const Box& extent = tetra.extent();
        if (extent.zmin() <= 0 && extent.zmax() >= 0) plotxy.writePolyhedron(coords, indices);
        if (extent.ymin() <= 0 && extent.ymax() >= 0) plotxz.writePolyhedron(coords, indices);
        if (extent.xmin() <= 0 && extent.xmax() >= 0) plotyz.writePolyhedron(coords, indices);
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
                const Tetra tetra = _grid->_tetrahedra[_mr];
                const array<Face, 4>& faces = tetra.faces();
                Position pos = r();
                Direction dir = k();

                int leavingFace = -1;
                double ds = DBL_MAX;

                // find entering face using a single Plücker product
                if (_enteringFace == -1) _enteringFace = tetra.findEnteringFace(pos, dir);

                // the translated Plücker moment in the local coordinate system
                Vec moment = Vec::cross(dir, pos - tetra.vertex(_enteringFace));

                // clockwise vertices around vertex 0
                std::array<int, 3> cv = clockwiseVertices(_enteringFace);

                // determine orientations for use in the decision tree
                double prod0 = Vec::dot(moment, tetra.edge(cv[0], _enteringFace));
                int clock0 = prod0 < 0;
                // if clockwise move clockwise else move cclockwise
                int i = clock0 ? 1 : 2;
                double prodi = Vec::dot(moment, tetra.edge(cv[i], _enteringFace));
                int cclocki = prodi >= 0;

                // use plane intersection algorithm if Plücker products are ambiguous
                // this is actually more strict than the algorithm described by Maria (2017)
                // but these edge cases are incredibly rare and can cause issues
                if (prod0 == 0. || prodi == 0.)
                {
                    for (int face : cv)
                    {
                        const Vec& n = faces[face]._normal;
                        double ndotk = Vec::dot(n, dir);
                        if (ndotk > 0)
                        {
                            const Vec& v = tetra.vertex(_enteringFace);
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
                    const Vec& n = faces[leavingFace]._normal;
                    const Vec& v = tetra.vertex(_enteringFace);
                    double ndotk = Vec::dot(n, dir);
                    ds = Vec::dot(n, v - pos) / ndotk;
                }

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
                    _mr = faces[leavingFace]._ntetra;

                    if (_mr < 0)
                    {
                        setState(State::Outside);
                        return false;
                    }
                    else
                    {
                        _enteringFace = faces[leavingFace]._nface;
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
