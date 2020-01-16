/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "DustTemperaturePerCellProbe.hpp"
#include "Configuration.hpp"
#include "MediumSystem.hpp"
#include "SpatialGrid.hpp"
#include "StringUtils.hpp"
#include "TextOutFile.hpp"
#include "Units.hpp"
#include "WavelengthGrid.hpp"

////////////////////////////////////////////////////////////////////

void DustTemperaturePerCellProbe::probeRun()
{
    if (find<Configuration>()->hasPanRadiationField() && find<MediumSystem>()->hasDust())
    {
        // locate the medium system
        auto ms = find<MediumSystem>();
        auto units = find<Units>();

        // create a text file
        TextOutFile file(this, itemName() + "_T", "dust temperature per cell");

        // write the header
        file.writeLine("# Indicative dust temperature per spatial cell");
        file.addColumn("spatial cell index", "", 'd');
        file.addColumn("indicative dust temperature", units->utemperature(), 'g');

        // write a line for each cell
        int numCells = ms->numCells();
        for (int m = 0; m != numCells; ++m)
        {
            file.writeRow(m, units->otemperature(ms->indicativeDustTemperature(m)));
        }
    }
}

////////////////////////////////////////////////////////////////////
