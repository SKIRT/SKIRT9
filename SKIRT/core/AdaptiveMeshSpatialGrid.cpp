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
#include "SpatialGridPlotFile.hpp"

////////////////////////////////////////////////////////////////////

void AdaptiveMeshSpatialGrid::setupSelfBefore()
{
    SpatialGrid::setupSelfBefore();

    // locate the adaptive mesh
    auto ms = find<MediumSystem>(false);
    auto interface = ms ? ms->interface<AdaptiveMeshInterface>() : nullptr;
    if (!interface) throw FATALERROR("Can't find an adaptive mesh in the medium system");
    _snapshot = interface->mesh();

    // tell it to construct neighbor information for its cells
    find<Log>()->info("Adding neighbor information to adaptive mesh...");
    _snapshot->addNeighbors();

    // if there is a single medium component, calculate the normalization factor imposed by it;
    // we need this to directly compute cell densities for the DensityInCellInterface
    if (ms->media().size() == 1)  _norm = ms->media()[0]->mass() / _snapshot->mass();
}

//////////////////////////////////////////////////////////////////////

int AdaptiveMeshSpatialGrid::dimension() const
{
    return 3;
}

//////////////////////////////////////////////////////////////////////

int AdaptiveMeshSpatialGrid::numCells() const
{
    return _snapshot->numEntities();
}

//////////////////////////////////////////////////////////////////////

Box AdaptiveMeshSpatialGrid::boundingBox() const
{
    return _snapshot->extent();
}

//////////////////////////////////////////////////////////////////////

double AdaptiveMeshSpatialGrid::volume(int m) const
{
    return _snapshot->volume(m);
}

//////////////////////////////////////////////////////////////////////

int AdaptiveMeshSpatialGrid::cellIndex(Position bfr) const
{
    return _snapshot->cellIndex(bfr);
}

//////////////////////////////////////////////////////////////////////

Position AdaptiveMeshSpatialGrid::centralPositionInCell(int m) const
{
    return _snapshot->position(m);
}

//////////////////////////////////////////////////////////////////////

Position AdaptiveMeshSpatialGrid::randomPositionInCell(int m) const
{
    return _snapshot->generatePosition(m);
}

//////////////////////////////////////////////////////////////////////

double AdaptiveMeshSpatialGrid::density(int /*h*/, int m) const
{
    return _norm * _snapshot->density(m);
}

//////////////////////////////////////////////////////////////////////

void AdaptiveMeshSpatialGrid::path(SpatialGridPath* path) const
{
    _snapshot->path(path);
}

//////////////////////////////////////////////////////////////////////

void AdaptiveMeshSpatialGrid::write_xy(SpatialGridPlotFile* outfile) const
{
    // Output the domain
    double xmin, ymin, zmin, xmax, ymax, zmax;
    boundingBox().extent(xmin, ymin, zmin, xmax, ymax, zmax);
    outfile->writeRectangle(xmin, ymin, xmax, ymax);

    // Output all leaf cells that intersect the coordinate plane
    double eps = 1e-8*(zmax-zmin);
    for (int m=0; m<numCells(); m++)
    {
        _snapshot->extent(m).extent(xmin, ymin, zmin, xmax, ymax, zmax);
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
    double eps = 1e-8*(ymax-ymin);
    for (int m=0; m<numCells(); m++)
    {
        _snapshot->extent(m).extent(xmin, ymin, zmin, xmax, ymax, zmax);
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
    double eps = 1e-8*(xmax-xmin);
    for (int m=0; m<numCells(); m++)
    {
        _snapshot->extent(m).extent(xmin, ymin, zmin, xmax, ymax, zmax);
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
    for (int m=0; m<numCells(); m++)
    {
        double xmin, ymin, zmin, xmax, ymax, zmax;
        _snapshot->extent(m).extent(xmin, ymin, zmin, xmax, ymax, zmax);
        outfile->writeCube(xmin, ymin, zmin, xmax, ymax, zmax);
    }
}

//////////////////////////////////////////////////////////////////////
