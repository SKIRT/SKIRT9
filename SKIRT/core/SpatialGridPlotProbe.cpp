/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "SpatialGridPlotProbe.hpp"
#include "MediumSystem.hpp"
#include "SpatialGrid.hpp"
#include "SpatialGridPlotFile.hpp"

////////////////////////////////////////////////////////////////////

void SpatialGridPlotProbe::probeSetup()
{
    // locate the grid (it is OK for the medium system to have no media components)
    auto ms = find<MediumSystem>(false);
    if (ms)
    {
        auto grid = ms->grid();
        int dimension = grid->dimension();

        // For the xy plane (always)
        {
            SpatialGridPlotFile outfile(this, itemName() + "_grid_xy");
            grid->write_xy(&outfile);
        }

        // For the xz plane (only if dimension is at least 2)
        if (dimension >= 2)
        {
            SpatialGridPlotFile outfile(this, itemName() + "_grid_xz");
            grid->write_xz(&outfile);
        }

        // For the yz plane (only if dimension is 3)
        if (dimension == 3)
        {
            SpatialGridPlotFile outfile(this, itemName() + "_grid_yz");
            grid->write_yz(&outfile);
        }

        // Full 3D coordinates (only if dimension is 3)
        if (dimension == 3)
        {
            SpatialGridPlotFile outfile(this, itemName() + "_grid_xyz");
            grid->write_xyz(&outfile);
        }
    }
}

////////////////////////////////////////////////////////////////////
