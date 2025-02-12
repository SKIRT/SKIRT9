/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "TetraMeshSpatialGrid.hpp"
#include "FatalError.hpp"
#include "Log.hpp"
#include "MediumSystem.hpp"
#include "NR.hpp"
#include "Random.hpp"
#include "SpatialGridPlotFile.hpp"
#include "StringUtils.hpp"
#include "TextInFile.hpp"
#include "tetgen.h"

//////////////////////////////////////////////////////////////////////

namespace
{
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

TetraMeshSpatialGrid::Tetra::Tetra(const vector<Vec>& vertices, const FourIndices& vertexIndices,
                                   const FourFaces& faces)
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
    static constexpr int etable[6][2] = {{3, 2}, {1, 3}, {2, 1}, {3, 0}, {0, 2}, {1, 0}};  // must be static
    int e = 0;
    // loop over all 6 edges because of rare cases where ray is inside edge and only 1 non-zero Plücker product
    for (int t1 = 0; t1 < 3; t1++)
    {
        for (int t2 = t1 + 1; t2 < 4; t2++)
        {
            Vec moment12 = Vec::cross(dir, pos - vertex(t1));
            double prod12 = Vec::dot(moment12, edge(t1, t2));
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

    if (Vec::dot(v - bfr, face._normal) < 0.) return false;

    return true;
}

//////////////////////////////////////////////////////////////////////

double TetraMeshSpatialGrid::Tetra::generateBarycentric(double& s, double& t, double& u) const
{
    if (s + t > 1.)  // cut'n fold the cube into a prism
    {
        s = 1. - s;
        t = 1. - t;
    }
    if (t + u > 1.)  // cut'n fold the prism into a tetrahedron
    {
        double tmp = u;
        u = 1. - s - t;
        t = 1. - tmp;
    }
    else if (s + t + u > 1.)
    {
        double tmp = u;
        u = s + t + u - 1.;
        s = 1. - t - tmp;
    }
    return 1. - u - t - s;
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
    return 1. / 6. * abs(Vec::dot(Vec::cross(edge(0, 1), edge(0, 2)), edge(0, 3)));
}

//////////////////////////////////////////////////////////////////////

double TetraMeshSpatialGrid::Tetra::diagonal() const
{
    double sum = 0.;
    for (int i = 0; i < 3; ++i)
    {
        for (int j = i + 1; j < 4; ++j)
        {
            sum += edge(i, j).norm2();
        }
    }
    return sqrt(sum / 6.);
}

//////////////////////////////////////////////////////////////////////

const TetraMeshSpatialGrid::FourFaces& TetraMeshSpatialGrid::Tetra::faces() const
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

void TetraMeshSpatialGrid::setupSelfBefore()
{
    BoxSpatialGrid::setupSelfBefore();

    _log = find<Log>();
    _eps = 1e-12 * extent().diagonal();

    // determine an appropriate set of samples and construct the Tetra mesh
    switch (_policy)
    {
        case Policy::Uniform:
        {
            auto random = find<Random>();
            _vertices.resize(_numSamples);
            for (int m = 0; m != _numSamples; ++m) _vertices[m] = random->position(extent());
            break;
        }
        case Policy::CentralPeak:
        {
            auto random = find<Random>();
            const int a = 1000;  // steepness of the peak; the central 1/a portion is NOT covered
            const double rscale = extent().rmax().norm();
            _vertices.resize(_numSamples);
            _vertices[0] = Vec();
            for (int m = 1; m != _numSamples;)  // skip first vertex so that it remains (0,0,0)
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
            _vertices = sampleMedia(media, weights, extent(), _numSamples);
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
            _vertices = sampleMedia(media, weights, extent(), _numSamples);
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
            _vertices = sampleMedia(media, weights, extent(), _numSamples);
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

    addCorners();
    removeInvalid();
    buildMesh();
    buildSearch();
}

////////////////////////////////////////////////////////////////////

void TetraMeshSpatialGrid::removeInvalid()
{
    _numVertices = _vertices.size();

    // remove vertices outside of the domain
    auto verticesEnd = std::remove_if(_vertices.begin(), _vertices.end(), [this](const Vec& vertex) {
        if (!contains(vertex)) return true;  // remove vertex
        return false;
    });
    _vertices.erase(verticesEnd, _vertices.end());

    // log removed vertices
    int numOutside = _numVertices - _vertices.size();
    if (numOutside) _log->info("removed " + StringUtils::toString(numOutside, 'd') + " vertices outside of the domain");
    _numVertices = _vertices.size();

    // sort vertices in order of increasing x coordinate to accelerate search for nearby sites
    std::sort(_vertices.begin(), _vertices.end(), [](Vec& v1, Vec& v2) { return v1.x() < v2.x(); });
    // mark vertices that lie too nearby another site
    std::vector<int> toRemove;
    for (int i = 0; i < _numVertices; ++i)
    {
        for (int j = i + 1; j < _numVertices && (_vertices[j].x() - _vertices[i].x() < _eps); ++j)
        {
            if ((_vertices[j] - _vertices[i]).norm2() < _eps * _eps)
            {
                toRemove.push_back(i);
                break;
            }
        }
    }
    // remove marked vertices in reverse order
    std::sort(toRemove.begin(), toRemove.end(), [](int i, int j) { return i > j; });
    for (int index : toRemove) _vertices.erase(_vertices.begin() + index);

    // log removed vertices
    int numNearby = toRemove.size();
    if (numNearby) _log->info("removed " + StringUtils::toString(numNearby, 'd') + " vertices too nearby to others");
    _numVertices = _vertices.size();
}

////////////////////////////////////////////////////////////////////

void TetraMeshSpatialGrid::addCorners()
{
    _numVertices = _vertices.size();
    // add the 8 corners of the domain to the list of vertices
    double xmin, ymin, zmin, xmax, ymax, zmax;
    extent(xmin, ymin, zmin, xmax, ymax, zmax);
    _vertices.reserve(_numVertices + 8);
    _vertices.emplace_back(xmin, ymin, zmin);
    _vertices.emplace_back(xmin, ymin, zmax);
    _vertices.emplace_back(xmin, ymax, zmin);
    _vertices.emplace_back(xmin, ymax, zmax);
    _vertices.emplace_back(xmax, ymin, zmin);
    _vertices.emplace_back(xmax, ymin, zmax);
    _vertices.emplace_back(xmax, ymax, zmin);
    _vertices.emplace_back(xmax, ymax, zmax);
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

    behavior.quiet = 1;  // -Q no console logging
    behavior.psc = 1;    // -s build Delaunay tetrahedralization
    // correct output options for out
    behavior.neighout = 2;   // -nn
    behavior.facesout = 1;   // -f
    behavior.zeroindex = 1;  // -z

    _log->info("Building Delaunay triangulation using " + StringUtils::toString(_numVertices, 'd')
               + " input vertices...");
    tetrahedralize(&behavior, &in, &out);
    _log->info("Built Delaunay triangulation");
}

////////////////////////////////////////////////////////////////////

void TetraMeshSpatialGrid::refineDelaunay(tetgenio& in, tetgenio& out)
{
    tetgenbehavior behavior;

    behavior.quiet = 1;  // -Q no console logging
    // tetgen refine options
    behavior.refine = 1;   // -r
    behavior.quality = 1;  // -q with default tetgen options for quality
    // correct output options for out
    behavior.neighout = 2;   // -nn
    behavior.facesout = 1;   // -f
    behavior.zeroindex = 1;  // -z

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
        // vertices
        FourIndices vertexIndices;
        for (int c = 0; c < 4; c++)
        {
            vertexIndices[c] = final.tetrahedronlist[4 * i + c];
        }

        // faces
        FourFaces faces;
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
            auto cv = clockwiseVertices(f);
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
    double avgVol = 1. / _numCells;
    double varVol = totalVol2 / _numCells / (V * V) - avgVol * avgVol;

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
    _log->info("Constructing search grid for " + std::to_string(_tetrahedra.size()) + " tetrahedra...");
    auto bounds = [this](int m) { return _tetrahedra[m].extent(); };
    auto intersects = [this](int m, const Box& box) { return box.intersects(_tetrahedra[m].extent()); };
    _search.loadEntities(_tetrahedra.size(), bounds, intersects);

    int nb = _search.numBlocks();
    _log->info("  Number of blocks in grid: " + std::to_string(nb * nb * nb) + " (" + std::to_string(nb) + "^3)");
    _log->info("  Smallest number of tetrahedra per block: " + std::to_string(_search.minEntitiesPerBlock()));
    _log->info("  Largest  number of tetrahedra per block: " + std::to_string(_search.maxEntitiesPerBlock()));
    _log->info("  Average  number of tetrahedra per block: "
               + StringUtils::toString(_search.avgEntitiesPerBlock(), 'f', 1));
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
    for (int m : _search.entitiesFor(bfr))
    {
        if (_tetrahedra[m].contains(bfr)) return m;
    }
    return -1;
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
            auto faceIndices = clockwiseVertices(v);
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
        if (numDone % 2000 == 0) _log->infoIfElapsed("Computed tetrahedra: ", 2000);
    }
}

//////////////////////////////////////////////////////////////////////

class TetraMeshSpatialGrid::MySegmentGenerator : public PathSegmentGenerator
{
    const TetraMeshSpatialGrid* _grid{nullptr};
    int _mr{-1};
    int _enteringFace{-1};

public:
    MySegmentGenerator(const TetraMeshSpatialGrid* grid) : _grid(grid) {}

    bool next() override
    {
        if (state() == State::Unknown)
        {
            // try moving the photon packet inside the grid; if this is impossible, return an empty path
            if (!moveInside(_grid->extent(), _grid->_eps)) return false;

            // get the index of the cell containing the current position
            _mr = _grid->cellIndex(r());
            _enteringFace = -1;

            // very rare edge case where no cell is found at domain boundary
            if (_mr == -1) return false;

            // if the photon packet started outside the grid, return the corresponding nonzero-length segment;
            // otherwise fall through to determine the first actual segment
            if (ds() > 0.) return true;
        }

        // inside convex hull
        if (state() == State::Inside)
        {
            // loop in case no exit point was found (which should happen only rarely)
            while (true)
            {
                const Tetra tetra = _grid->_tetrahedra[_mr];
                const FourFaces& faces = tetra.faces();
                Position pos = r();
                Direction dir = k();

                int leavingFace = -1;
                double ds = DBL_MAX;

                // find entering face using a single Plücker product
                if (_enteringFace == -1) _enteringFace = tetra.findEnteringFace(pos, dir);

                // the translated Plücker moment in the local coordinate system
                Vec moment = Vec::cross(dir, pos - tetra.vertex(_enteringFace));

                // clockwise vertices around entering face
                auto cv = clockwiseVertices(_enteringFace);

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
                    static constexpr int dtable[2][2] = {{1, 0}, {0, 2}};  // must be static
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
