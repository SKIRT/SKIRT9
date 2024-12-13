/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "TetraMeshSpatialGrid.hpp"
#include "FatalError.hpp"
#include "MediumSystem.hpp"
#include "NR.hpp"
#include "Random.hpp"
#include "SpatialGridPlotFile.hpp"
#include "StringUtils.hpp"
#include "TextInFile.hpp"
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

    // get the vertices around a face in clockwise order viewed from opposite vertex
    inline std::array<int, 3> clockwiseVertices(int face)
    {
        std::array<int, 3> cv = {(face + 3) % 4, (face + 2) % 4, (face + 1) % 4};
        // if face is even we should swap two edges
        if (face % 2 == 0) std::swap(cv[0], cv[2]);
        return cv;
    }
}

//////////////////////////////////////////////////////////////////////

TetraMeshSpatialGrid::Tetra::Tetra(const vector<Vec>& vertices, const array<int, 4>& vertexIndices,
                                   const array<Face, 4>& faces)
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

//////////////////////////////////////////////////////////////////////

int TetraMeshSpatialGrid::Tetra::findEnteringFace(const Vec& pos, const Direction& dir) const
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

//////////////////////////////////////////////////////////////////////

bool TetraMeshSpatialGrid::Tetra::contains(const Position& bfr) const
{
    if (!_extent.contains(bfr)) return false;

    // Use the same vertex for the first 3 faces
    Vec v = vertex(3);  // vertex 3 is on faces 0,1,2

    for (int f = 0; f <= 2; f++)
    {
        const Face& face = _faces[f];

        // if bfr->v is opposite direction as the outward pointing normal, the point is outside
        if (Vec::dot(v - bfr, face._normal) < 0) return false;
    }

    // Check face 3
    const Face& face = _faces[3];
    v = vertex(0);  // any vertex that is not vertex 3

    if (Vec::dot(v - bfr, face._normal) < 0) return false;

    return true;
}

//////////////////////////////////////////////////////////////////////

double TetraMeshSpatialGrid::Tetra::generateBarycentric(double& s, double& t, double& u) const
{
    if (s + t > 1.0)  // cut'n fold the cube into a prism
    {
        s = 1.0 - s;
        t = 1.0 - t;
    }
    if (t + u > 1.0)  // cut'n fold the prism into a tetrahedron
    {
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

//////////////////////////////////////////////////////////////////////

Position TetraMeshSpatialGrid::Tetra::generatePosition(Random* random) const
{
    double s = random->uniform();
    double t = random->uniform();
    double u = random->uniform();

    double r = generateBarycentric(s, t, u);

    return Position(r * vertex(0) + u * vertex(1) + t * vertex(2) + s * vertex(3));
}

//////////////////////////////////////////////////////////////////////

Vec TetraMeshSpatialGrid::Tetra::vertex(int t) const
{
    return _vertices[_vertexIndices[t]];
}

//////////////////////////////////////////////////////////////////////

Vec TetraMeshSpatialGrid::Tetra::edge(int t1, int t2) const
{
    return vertex(t2) - vertex(t1);
}

//////////////////////////////////////////////////////////////////////

double TetraMeshSpatialGrid::Tetra::volume() const
{
    return 1 / 6. * abs(Vec::dot(Vec::cross(edge(0, 1), edge(0, 2)), edge(0, 3)));
}

//////////////////////////////////////////////////////////////////////

double TetraMeshSpatialGrid::Tetra::diagonal() const
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

//////////////////////////////////////////////////////////////////////

const array<TetraMeshSpatialGrid::Face, 4>& TetraMeshSpatialGrid::Tetra::faces() const
{
    return _faces;
}

//////////////////////////////////////////////////////////////////////

const Vec& TetraMeshSpatialGrid::Tetra::centroid() const
{
    return _centroid;
}

//////////////////////////////////////////////////////////////////////

const Box& TetraMeshSpatialGrid::Tetra::extent() const
{
    return _extent;
}

//////////////////////////////////////////////////////////////////////

class TetraMeshSpatialGrid::BlockGrid
{
    const vector<Tetra>& _tetrahedra;  // reference to list of all tetrahedra
    int _gridsize;                     // number of grid blocks in each spatial direction
    Array _xgrid, _ygrid, _zgrid;      // the m+1 grid separation points for each spatial direction
    vector<vector<int>> _listv;        // the m*m*m lists of indices for cells overlapping each grid block
    int _pmin, _pmax;                  // statistics of the minimum and maximum number of cells per block

public:
    // The constructor creates a cuboidal grid of the specified number of grid blocks in each
    // spatial direction, and for each of the grid blocks it builds a list of all blocks
    // overlapping the grid block. In an attempt to distribute the blocks evenly over the
    // grid blocks, the sizes of the grid blocks in each spatial direction are chosen so that
    // the block centers are evenly distributed over the grid blocks.
    BlockGrid(const vector<Tetra>& tetrahedra, Box extent, int gridsize) : _tetrahedra(tetrahedra), _gridsize(gridsize)
    {
        // build the grids in each spatial direction
        makegrid(0, gridsize, _xgrid, extent.xmin(), extent.xmax());
        makegrid(1, gridsize, _ygrid, extent.ymin(), extent.ymax());
        makegrid(2, gridsize, _zgrid, extent.zmin(), extent.zmax());

        // make room for p*p*p grid blocks
        _listv.resize(gridsize * gridsize * gridsize);

        // add each block to the list for every grid block that it overlaps
        int n = _tetrahedra.size();
        for (int m = 0; m != n; ++m)
        {
            Box boundingBox = _tetrahedra[m].extent();

            // find indices for first and last grid block overlapped by block, in each spatial direction
            int i1 = NR::locateClip(_xgrid, boundingBox.xmin());
            int i2 = NR::locateClip(_xgrid, boundingBox.xmax());
            int j1 = NR::locateClip(_ygrid, boundingBox.ymin());
            int j2 = NR::locateClip(_ygrid, boundingBox.ymax());
            int k1 = NR::locateClip(_zgrid, boundingBox.zmin());
            int k2 = NR::locateClip(_zgrid, boundingBox.zmax());

            // add the block to all grid blocks in that 3D range
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
        for (int index = 0; index < gridsize * gridsize * gridsize; index++)
        {
            int size = _listv[index].size();
            _pmin = min(_pmin, size);
            _pmax = max(_pmax, size);
        }
    }

    // This function determines the grid separation points along a specified axis (x, y, or z)
    // to ensure that the cells are evenly distributed across the grid blocks. It does this by
    // binning the tetrahedra centers at a high resolution and then calculating the cumulative
    // distribution to set the grid separation points.
    void makegrid(int axis, int gridsize, Array& grid, double cmin, double cmax)
    {
        int n = _tetrahedra.size();

        // determine the block distribution by binning at a decent resolution
        int nbins = gridsize * 100;
        double binwidth = (cmax - cmin) / nbins;
        vector<int> bins(nbins);
        for (const Tetra& tetra : _tetrahedra)
        {
            double center;
            switch (axis)
            {
                case 0: center = tetra.centroid().x(); break;
                case 1: center = tetra.centroid().y(); break;
                case 2: center = tetra.centroid().z(); break;
            }
            bins[static_cast<int>((center - cmin) / binwidth)] += 1;
        }

        // determine grid separation points based on the cumulative distribution
        grid.resize(gridsize + 1);
        grid[0] = -std::numeric_limits<double>::infinity();
        int perblock = n / gridsize;  // target number of particles per block
        int cumul = 0;                // cumulative number of particles in processed bins
        int gridindex = 1;            // index of the next grid separation point to be filled
        for (int binindex = 0; binindex < nbins; binindex++)
        {
            cumul += bins[binindex];
            if (cumul > perblock * gridindex)
            {
                grid[gridindex] = cmin + (binindex + 1) * binwidth;
                gridindex += 1;
                if (gridindex >= gridsize) break;
            }
        }
        grid[gridsize] = std::numeric_limits<double>::infinity();
    }

    // This function returns the smallest number of blocks overlapping a single grid block.
    int minCellRefsPerBlock() const { return _pmin; }

    // This function returns the largest number of blocks overlapping a single grid block.
    int maxCellRefsPerBlock() const { return _pmax; }

    // This function returns the index (in the list originally passed to the constructor)
    // of the first block in the list that overlaps the specified position,
    // or -1 if none of the blocks in the list overlap the specified position.
    int cellIndexFor(Position r) const
    {
        // locate the grid block containing the specified position
        int i = NR::locateClip(_xgrid, r.x());
        int j = NR::locateClip(_ygrid, r.y());
        int k = NR::locateClip(_zgrid, r.z());

        // search the list of blocks for that grid block
        for (int m : _listv[index(_gridsize, i, j, k)])
        {
            if (_tetrahedra[m].contains(r)) return m;
        }
        return -1;
    }
};

////////////////////////////////////////////////////////////////////

TetraMeshSpatialGrid::~TetraMeshSpatialGrid()
{
    delete _blocks;
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
        case Policy::File:
        {
            // read the input file
            TextInFile in(this, _filename, "tetrahedral vertices");
            in.addColumn("position x", "length", "pc");
            in.addColumn("position y", "length", "pc");
            in.addColumn("position z", "length", "pc");
            Array coords;
            while (in.readRow(coords)) _vertices.emplace_back(coords[0], coords[1], coords[2]);
            in.close();
        }
    }

    removeOutside();
    if (_numVertices < 4) return;  // abort if there are not enough vertices
    buildMesh();
    buildSearch();
}

////////////////////////////////////////////////////////////////////

void TetraMeshSpatialGrid::removeOutside()
{
    _numVertices = _vertices.size();

    // remove vertices outside of the domain
    auto sitesEnd = std::remove_if(_vertices.begin(), _vertices.end(), [this](const Vec& vertex) {
        if (!contains(vertex)) return true;  // remove vertex
        return false;
    });
    _vertices.erase(sitesEnd, _vertices.end());

    // log removed vertices
    int numOutside = _numVertices - _vertices.size();
    if (numOutside) _log->info("removed " + StringUtils::toString(numOutside, 'd') + " vertices outside of the domain");

    _numVertices = _vertices.size();
}

////////////////////////////////////////////////////////////////////

void TetraMeshSpatialGrid::buildMesh()
{
    tetgenio delaunay;
    buildDelaunay(delaunay);

    if (refine())
    {
        tetgenio refined;
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
    tetgenio in;
    tetgenbehavior behavior;

    in.numberofpoints = _vertices.size();
    in.pointlist = new REAL[in.numberofpoints * 3];
    for (int i = 0; i < in.numberofpoints; i++)
    {
        in.pointlist[i * 3 + 0] = _vertices[i].x();
        in.pointlist[i * 3 + 1] = _vertices[i].y();
        in.pointlist[i * 3 + 2] = _vertices[i].z();
    }

    behavior.psc = 1;  // -s build Delaunay tetrahedralization
    // correct output options for out
    behavior.neighout = 2;   // -nn
    behavior.facesout = 1;   // -f
    behavior.zeroindex = 1;  // -z

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

    behavior.quiet = 1;  // -Q

    _log->info("Refining triangulation...");
    tetrahedralize(&behavior, &in, &out);
    _log->info("Built refined triangulation");
}

////////////////////////////////////////////////////////////////////

void TetraMeshSpatialGrid::storeTetrahedra(const tetgenio& final, bool storeVertices)
{
    _numCells = final.numberoftetrahedra;

    // replace old vertices
    if (storeVertices)
    {
        _numVertices = final.numberofpoints;

        _vertices.resize(_numVertices);
        for (int i = 0; i < _numVertices; i++)
        {
            double x = final.pointlist[3 * i + 0];
            double y = final.pointlist[3 * i + 1];
            double z = final.pointlist[3 * i + 2];

            _vertices[i] = Vec(x, y, z);
        }
    }

    _tetrahedra.reserve(_numCells);  // no default constructor for Tetra
    for (int i = 0; i < _numCells; i++)
    {
        std::array<int, 4> vertexIndices;
        std::array<Face, 4> faces;

        // vertices
        for (int c = 0; c < 4; c++)
        {
            vertexIndices[c] = final.tetrahedronlist[4 * i + c];
        }

        // faces
        for (int f = 0; f < 4; f++)
        {
            // -1 if no neighbor
            int ntetra = final.neighborlist[4 * i + f];

            // find which face is shared with neighbor
            int nface = -1;
            if (ntetra != -1)
            {
                for (int fn = 0; fn < 4; fn++)
                {
                    if (final.neighborlist[4 * ntetra + fn] == i)
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

    // log statistics
    _log->info("Done computing tetrahedralization");
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
    _log->info("Constructing intermediate " + size + "x" + size + "x" + size + " grid for tetrahedra...");
    _blocks = new BlockGrid(_tetrahedra, extent(), gridsize);
    _log->info("  Smallest number of tetrahedra per grid block: " + std::to_string(_blocks->minCellRefsPerBlock()));
    _log->info("  Largest  number of tetrahedra per grid block: " + std::to_string(_blocks->maxCellRefsPerBlock()));
    _log->info("  Average  number of tetrahedra per grid block: "
               + StringUtils::toString(_numCells / double(gridsize * gridsize * gridsize), 'f', 1));
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
    return _blocks->cellIndexFor(bfr);
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

    // for each tetrahedron, compute its edges
    _log->info("Writing plot files for tetrahedral grid with " + std::to_string(_numCells) + " tetrahedra");
    _log->infoSetElapsed(_numCells);
    int numDone = 0;
    for (int i = 0; i < _numCells; i++)
    {
        const Tetra& tetra = _tetrahedra[i];
        vector<double> coords(12);
        vector<int> indices(16);

        // write each face as a polygon
        for (int v = 0; v < 4; v++)
        {
            Vec vertex = tetra.vertex(v);
            coords[3 * v + 0] = vertex.x();
            coords[3 * v + 1] = vertex.y();
            coords[3 * v + 2] = vertex.z();

            // get vertices of opposite face
            std::array<int, 3> faceIndices = clockwiseVertices(v);
            indices[4 * v + 0] = 3;  // amount of vertices per face
            indices[4 * v + 1] = faceIndices[0];
            indices[4 * v + 2] = faceIndices[1];
            indices[4 * v + 3] = faceIndices[2];
        }

        const Box& extent = tetra.extent();
        if (extent.zmin() <= 0 && extent.zmax() >= 0) plotxy.writePolyhedron(coords, indices);
        if (extent.ymin() <= 0 && extent.ymax() >= 0) plotxz.writePolyhedron(coords, indices);
        if (extent.xmin() <= 0 && extent.xmax() >= 0) plotyz.writePolyhedron(coords, indices);
        if (i <= 1000) plotxyz.writePolyhedron(coords, indices);

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
            // moveInside does not move the photon packet inside the convex hull so it is currently disabled
            // if (!moveInside(_grid->extent(), _grid->_eps)) return false;

            // get the index of the cell containing the current position
            _mr = _grid->cellIndex(r());
            _enteringFace = -1;

            // position is outside of convex hull
            if (_mr < 0) return false;
        }

        // inside convex hull
        if (state() == State::Inside)
        {
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

                // clockwise vertices around entering face
                std::array<int, 3> cv = clockwiseVertices(_enteringFace);

                // determine orientations for use in the decision tree
                double prod0 = Vec::dot(moment, tetra.edge(cv[0], _enteringFace));
                int clock0 = prod0 < 0;
                // if clockwise move clockwise else move counterclockwise
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
                    // plane intersection to leaving face
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
