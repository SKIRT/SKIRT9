/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "VoronoiMeshSpatialGrid.hpp"
#include "FatalError.hpp"
#include "Log.hpp"
#include "MediumSystem.hpp"
#include "Random.hpp"
#include "SiteListInterface.hpp"
#include "SpatialGridPlotFile.hpp"
#include "VoronoiMeshInterface.hpp"
#include "VoronoiMeshSnapshot.hpp"
#include "container.hh"

//////////////////////////////////////////////////////////////////////

VoronoiMeshSpatialGrid::~VoronoiMeshSpatialGrid()
{
    if (_policy != Policy::ImportedMesh) delete _mesh;
}

//////////////////////////////////////////////////////////////////////

void VoronoiMeshSpatialGrid::setupSelfBefore()
{
    BoxSpatialGrid::setupSelfBefore();

    // determine an appropriate set of sites and construct the Voronoi mesh
    switch (_policy)
    {
    case Policy::Uniform:
        {
            auto random = find<Random>();
            vector<Vec> rv(_numSites);
            for (int m=0; m!=_numSites; ++m) rv[m] = random->position(extent());
            _mesh = new VoronoiMeshSnapshot(this, extent(), rv);
            break;
        }
    case Policy::CentralPeak:
        {
            auto random = find<Random>();
            const int a = 1000;                     // steepness of the peak; the central 1/a portion is NOT covered
            const double rscale = extent().rmax().norm();
            vector<Vec> rv(_numSites);
            for (int m=1; m!=_numSites; ++m)        // skip first particle so that it remains (0,0,0)
            {
                while (true)
                {
                    double r = rscale * pow(1./a, random->uniform());   // random distribution according to 1/x
                    Direction k = random->direction();
                    Position p = Position(r,k);
                    if (extent().contains(p))       // discard any points outside of the domain
                    {
                        rv[m] = p;
                        break;
                    }
                }
            }
            _mesh = new VoronoiMeshSnapshot(this, extent(), rv);
            break;
        }
    case Policy::MediumDensity:
        {
            auto medium = find<MediumSystem>()->find<Medium>();
            vector<Vec> rv(_numSites);
            for (int m=0; m!=_numSites; ++m)
            {
                while (true)
                {
                    Position p = medium->generatePosition();
                    if (extent().contains(p))       // discard any points outside of the domain
                    {
                        rv[m] = p;
                        break;
                    }
                }
            }
            _mesh = new VoronoiMeshSnapshot(this, extent(), rv);
            break;
        }
    case Policy::File:
        {
            _mesh = new VoronoiMeshSnapshot(this, extent(), _filename);
            break;
        }
    case Policy::ImportedSites:
        {
            auto sli = find<MediumSystem>()->interface<SiteListInterface>(2);
            _mesh = new VoronoiMeshSnapshot(this, extent(), sli);
            break;
        }
    case Policy::ImportedMesh:
        {
            auto ms = find<MediumSystem>(false);
            _mesh = ms->interface<VoronoiMeshInterface>(2)->voronoiMesh();

            // if there is a single medium component, calculate the normalization factor imposed by it;
            // we need this to directly compute cell densities for the DensityInCellInterface
            if (ms->media().size() == 1)  _norm = ms->media()[0]->number() / _mesh->mass();
            break;
        }
    }
}

//////////////////////////////////////////////////////////////////////

int VoronoiMeshSpatialGrid::numCells() const
{
    return _mesh->numEntities();
}

//////////////////////////////////////////////////////////////////////

double VoronoiMeshSpatialGrid::volume(int m) const
{
    return _mesh->volume(m);
}

//////////////////////////////////////////////////////////////////////

int VoronoiMeshSpatialGrid::cellIndex(Position bfr) const
{
    return _mesh->cellIndex(bfr);
}

//////////////////////////////////////////////////////////////////////

Position VoronoiMeshSpatialGrid::centralPositionInCell(int m) const
{
    return _mesh->centroidPosition(m);
}

//////////////////////////////////////////////////////////////////////

Position VoronoiMeshSpatialGrid::randomPositionInCell(int m) const
{
    return _mesh->generatePosition(m);
}

//////////////////////////////////////////////////////////////////////

void VoronoiMeshSpatialGrid::path(SpatialGridPath* path) const
{
    _mesh->path(path);
}

//////////////////////////////////////////////////////////////////////

void VoronoiMeshSpatialGrid::performWriteGrid() const
{
    // create the plot files
    SpatialGridPlotFile plotxy(this, "ds_gridxy");
    SpatialGridPlotFile plotxz(this, "ds_gridxz");
    SpatialGridPlotFile plotyz(this, "ds_gridyz");
    SpatialGridPlotFile plotxyz(this, "ds_gridxyz");

    // load all sites in a Voro container
    int numCells = _mesh->numEntities();
    int nb = max(3, min(1000, static_cast<int>(pow(numCells/5.,1./3.)) ));
    voro::container con(xmin(), xmax(), ymin(), ymax(), zmin(), zmax(), nb, nb, nb, false,false,false, 8);
    for (int m=0; m!=numCells; ++m)
    {
        Vec r = _mesh->position(m);
        con.put(m, r.x(),r.y(),r.z());
    }

    // loop over all Voro cells
    voro::c_loop_all loop(con);
    if (loop.start()) do
    {
        // Compute the cell
        voro::voronoicell fullcell;
        con.compute_cell(fullcell, loop);

        // Get the edges of the cell
        double x,y,z;
        loop.pos(x,y,z);
        vector<double> coords;
        fullcell.vertices(x,y,z, coords);
        vector<int> indices;
        fullcell.face_vertices(indices);

        // Write the edges of the cell to the plot files
        Box bounds = _mesh->extent(loop.pid());
        if (bounds.zmin()<=0 && bounds.zmax()>=0) plotxy.writePolyhedron(coords, indices);
        if (bounds.ymin()<=0 && bounds.ymax()>=0) plotxz.writePolyhedron(coords, indices);
        if (bounds.xmin()<=0 && bounds.xmax()>=0) plotyz.writePolyhedron(coords, indices);
        if (loop.pid() <= 1000) plotxyz.writePolyhedron(coords, indices);
    }
    while (loop.inc());
}

//////////////////////////////////////////////////////////////////////

double VoronoiMeshSpatialGrid::numberDensity(int /*h*/, int m) const
{
    return _norm * _mesh->density(m);
}

//////////////////////////////////////////////////////////////////////

bool VoronoiMeshSpatialGrid::offersInterface(const std::type_info& interfaceTypeInfo) const
{
    if (interfaceTypeInfo == typeid(DensityInCellInterface))
    {
        if (_policy != Policy::ImportedMesh) return false;
        auto ms = find<MediumSystem>(false);
        return ms && ms->media().size() == 1;
    }
    return SpatialGrid::offersInterface(interfaceTypeInfo);
}

//////////////////////////////////////////////////////////////////////
