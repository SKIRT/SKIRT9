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
#include "VoronoiMeshInterface.hpp"
#include "VoronoiMeshSnapshot.hpp"

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
            _mesh = new VoronoiMeshSnapshot(this, extent(), rv, _relaxSites);
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
            _mesh = new VoronoiMeshSnapshot(this, extent(), rv, _relaxSites);
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
            _mesh = new VoronoiMeshSnapshot(this, extent(), rv, _relaxSites);
            break;
        }
    case Policy::File:
        {
            _mesh = new VoronoiMeshSnapshot(this, extent(), _filename, _relaxSites);
            break;
        }
    case Policy::ImportedSites:
        {
            auto sli = find<MediumSystem>()->interface<SiteListInterface>(2);
            _mesh = new VoronoiMeshSnapshot(this, extent(), sli, _relaxSites);
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

void VoronoiMeshSpatialGrid::writeGridPlotFiles(const SimulationItem* probe) const
{
    _mesh->writeGridPlotFiles(probe);
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
