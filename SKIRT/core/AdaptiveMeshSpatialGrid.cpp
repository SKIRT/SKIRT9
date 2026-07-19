/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "AdaptiveMeshSpatialGrid.hpp"
#include "AdaptiveMeshInterface.hpp"
#include "AdaptiveMeshSnapshot.hpp"
#include "FatalError.hpp"
#include "Log.hpp"
#include "MediumSystem.hpp"
#include "PathSegmentGenerator.hpp"
#include "SpatialGridPlotFile.hpp"

////////////////////////////////////////////////////////////////////

void AdaptiveMeshSpatialGrid::setupSelfBefore()
{
    SpatialGrid::setupSelfBefore();

    // locate the adaptive mesh, scanning medium components in configuration order
    auto ms = find<MediumSystem>(false);
    if (!ms) throw FATALERROR("There is no medium system in which to locate an adaptive mesh");
    _mesh = ms->interface<AdaptiveMeshInterface>(2)->adaptiveMesh();

    // tell it to construct neighbor information for its cells
    find<Log>()->info("Adding neighbor information to adaptive mesh...");
    _mesh->addNeighbors();

    // if there is a single medium component, calculate the normalization factor imposed by it;
    // we need this to directly compute cell densities for the DensityInCellInterface
    if (ms->media().size() == 1) _norm = ms->media()[0]->number() / _mesh->mass();
}

//////////////////////////////////////////////////////////////////////

int AdaptiveMeshSpatialGrid::dimension() const
{
    return 3;
}

//////////////////////////////////////////////////////////////////////

int AdaptiveMeshSpatialGrid::numCells() const
{
    return _mesh->numEntities();
}

//////////////////////////////////////////////////////////////////////

Box AdaptiveMeshSpatialGrid::boundingBox() const
{
    return _mesh->extent();
}

//////////////////////////////////////////////////////////////////////

double AdaptiveMeshSpatialGrid::volume(int m) const
{
    return _mesh->volume(m);
}

//////////////////////////////////////////////////////////////////////

double AdaptiveMeshSpatialGrid::diagonal(int m) const
{
    return _mesh->diagonal(m);
}

//////////////////////////////////////////////////////////////////////

int AdaptiveMeshSpatialGrid::cellIndex(Position bfr) const
{
    return _mesh->cellIndex(bfr);
}

//////////////////////////////////////////////////////////////////////

Position AdaptiveMeshSpatialGrid::centralPositionInCell(int m) const
{
    return _mesh->position(m);
}

//////////////////////////////////////////////////////////////////////

Position AdaptiveMeshSpatialGrid::randomPositionInCell(int m) const
{
    return _mesh->generatePosition(m);
}

//////////////////////////////////////////////////////////////////////

std::unique_ptr<PathSegmentGenerator> AdaptiveMeshSpatialGrid::createPathSegmentGenerator() const
{
    return _mesh->createPathSegmentGenerator();
}

//////////////////////////////////////////////////////////////////////

void AdaptiveMeshSpatialGrid::write_xy(SpatialGridPlotFile* outfile) const
{
    // Output the domain
    double xmin, ymin, zmin, xmax, ymax, zmax;
    boundingBox().extent(xmin, ymin, zmin, xmax, ymax, zmax);
    outfile->writeRectangle(xmin, ymin, xmax, ymax);

    // Output all leaf cells that intersect the coordinate plane
    double eps = 1e-8 * (zmax - zmin);
    for (int m = 0; m < numCells(); m++)
    {
        _mesh->extent(m).extent(xmin, ymin, zmin, xmax, ymax, zmax);
        if (zmin < eps && zmax > -eps)
        {
            outfile->writeRectangle(xmin, ymin, xmax, ymax);
        }
    }
}

//////////////////////////////////////////////////////////////////////

void AdaptiveMeshSpatialGrid::write_xz(SpatialGridPlotFile* outfile) const
{
    // Output the domain
    double xmin, ymin, zmin, xmax, ymax, zmax;
    boundingBox().extent(xmin, ymin, zmin, xmax, ymax, zmax);
    outfile->writeRectangle(xmin, zmin, xmax, zmax);

    // Output all leaf cells that intersect the coordinate plane
    double eps = 1e-8 * (ymax - ymin);
    for (int m = 0; m < numCells(); m++)
    {
        _mesh->extent(m).extent(xmin, ymin, zmin, xmax, ymax, zmax);
        if (ymin < eps && ymax > -eps)
        {
            outfile->writeRectangle(xmin, zmin, xmax, zmax);
        }
    }
}

//////////////////////////////////////////////////////////////////////

void AdaptiveMeshSpatialGrid::write_yz(SpatialGridPlotFile* outfile) const
{
    // Output the domain
    double xmin, ymin, zmin, xmax, ymax, zmax;
    boundingBox().extent(xmin, ymin, zmin, xmax, ymax, zmax);
    outfile->writeRectangle(ymin, zmin, ymax, zmax);

    // Output all leaf cells that intersect the coordinate plane
    double eps = 1e-8 * (xmax - xmin);
    for (int m = 0; m < numCells(); m++)
    {
        _mesh->extent(m).extent(xmin, ymin, zmin, xmax, ymax, zmax);
        if (xmin < eps && xmax > -eps)
        {
            outfile->writeRectangle(ymin, zmin, ymax, zmax);
        }
    }
}

//////////////////////////////////////////////////////////////////////

void AdaptiveMeshSpatialGrid::write_xyz(SpatialGridPlotFile* outfile) const
{
    // Output all leaf cells (should we limit this somehow?)
    for (int m = 0; m < numCells(); m++)
    {
        double xmin, ymin, zmin, xmax, ymax, zmax;
        _mesh->extent(m).extent(xmin, ymin, zmin, xmax, ymax, zmax);
        outfile->writeCube(xmin, ymin, zmin, xmax, ymax, zmax);
    }
}

//////////////////////////////////////////////////////////////////////

double AdaptiveMeshSpatialGrid::numberDensity(int /*h*/, int m) const
{
    return _norm * _mesh->density(m);
}

//////////////////////////////////////////////////////////////////////

bool AdaptiveMeshSpatialGrid::offersInterface(const std::type_info& interfaceTypeInfo) const
{
    if (interfaceTypeInfo == typeid(DensityInCellInterface))
    {
        auto ms = find<MediumSystem>(false);
        return ms && ms->media().size() == 1;
    }
    return SpatialGrid::offersInterface(interfaceTypeInfo);
}

//////////////////////////////////////////////////////////////////////
