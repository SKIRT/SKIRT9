/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "SpatialGridSourceDensityProbe.hpp"
#include "GeometricSource.hpp"
#include "MediumSystem.hpp"
#include "SourceSystem.hpp"
#include "TextOutFile.hpp"
#include "Units.hpp"

////////////////////////////////////////////////////////////////////

void SpatialGridSourceDensityProbe::probe()
{
    // locate the grid (it is OK for the medium system to have no media components)
    auto ms = find<MediumSystem>(false);
    auto grid = ms ? ms->grid() : nullptr;

    // locate the geometric sources
    vector<int> hv;           // source component indices
    vector<Geometry*> geomv;  // corresponding geometries
    int h = 0;
    for (auto source : find<SourceSystem>()->sources())
    {
        auto geomsource = dynamic_cast<GeometricSource*>(source);
        if (geomsource)
        {
            hv.push_back(h);
            geomv.push_back(geomsource->geometry());
        }
        h++;
    }

    if (grid && hv.size())
    {
        auto units = find<Units>();

        // create a text file and add the columns
        TextOutFile out(this, itemName() + "_sourcedens", "gridded primary source densities");
        out.addColumn("spatial cell index", "", 'd');
        for (int h : hv)
            out.addColumn("normalized density for source " + std::to_string(h + 1), "1/" + units->uvolume());

        // write a line for each cell
        int numCells = grid->numCells();
        for (int m = 0; m != numCells; ++m)
        {
            Position p = grid->centralPositionInCell(m);
            vector<double> row;
            row.push_back(static_cast<double>(m));
            for (auto geom : geomv) row.push_back(1. / units->ovolume(1. / geom->density(p)));
            out.writeRow(row);
        }
    }
}

////////////////////////////////////////////////////////////////////
