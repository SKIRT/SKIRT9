/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "TreeSpatialGridTopologyProbe.hpp"
#include "MediumSystem.hpp"
#include "TextOutFile.hpp"
#include "TreeSpatialGrid.hpp"

////////////////////////////////////////////////////////////////////

void TreeSpatialGridTopologyProbe::probe()
{
    // locate the grid (it is OK for the medium system to have no media components)
    auto ms = find<MediumSystem>(false);
    auto grid = ms ? ms->find<TreeSpatialGrid>(false) : nullptr;
    if (grid)
    {
        TextOutFile outfile(this, itemName() + "_treetop", "spatial tree grid topology");
        grid->writeTopology(&outfile);
    }
}

////////////////////////////////////////////////////////////////////
