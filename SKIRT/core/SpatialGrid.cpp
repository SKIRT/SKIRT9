/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "SpatialGrid.hpp"
#include "Random.hpp"
#include "SpatialGridPlotFile.hpp"

//////////////////////////////////////////////////////////////////////

void SpatialGrid::setupSelfBefore()
{
    SimulationItem::setupSelfBefore();

    _random = find<Random>();
}

//////////////////////////////////////////////////////////////////////

void SpatialGrid::writeGridPlotFiles(const SimulationItem* probe) const
{
    // For the xy plane (always)
    {
        SpatialGridPlotFile outfile(probe, probe->itemName() + "_grid_xy");
        write_xy(&outfile);
    }

    // For the xz plane (only if dimension is at least 2)
    if (dimension() >= 2)
    {
        SpatialGridPlotFile outfile(probe, probe->itemName() + "_grid_xz");
        write_xz(&outfile);
    }

    // For the yz plane (only if dimension is 3)
    if (dimension() == 3)
    {
        SpatialGridPlotFile outfile(probe, probe->itemName() + "_grid_yz");
        write_yz(&outfile);
    }

    // Full 3D coordinates (only if dimension is 3)
    if (dimension() == 3)
    {
        SpatialGridPlotFile outfile(probe, probe->itemName() + "_grid_xyz");
        write_xyz(&outfile);
    }
}

//////////////////////////////////////////////////////////////////////

void SpatialGrid::write_xy(SpatialGridPlotFile* /*outfile*/) const {}

//////////////////////////////////////////////////////////////////////

void SpatialGrid::write_xz(SpatialGridPlotFile* /*outfile*/) const {}

//////////////////////////////////////////////////////////////////////

void SpatialGrid::write_yz(SpatialGridPlotFile* /*outfile*/) const {}

//////////////////////////////////////////////////////////////////////

void SpatialGrid::write_xyz(SpatialGridPlotFile* /*outfile*/) const {}

//////////////////////////////////////////////////////////////////////
