/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "GasTemperaturePerCellProbe.hpp"
#include "Configuration.hpp"
#include "MediumSystem.hpp"
#include "TextOutFile.hpp"
#include "Units.hpp"

////////////////////////////////////////////////////////////////////

void GasTemperaturePerCellProbe::probeSetup()
{
    if (find<Configuration>()->hasLymanAlpha())
    {
        // locate the medium system
        auto ms = find<MediumSystem>();
        auto units = find<Units>();

        // create a text file
        TextOutFile file(this, itemName() + "_T", "gas temperature per cell");

        // write the header
        file.writeLine("# Gas temperature per spatial cell");
        file.addColumn("spatial cell index", "", 'd');
        file.addColumn("gas temperature", units->utemperature(), 'g');

        // write a line for each cell
        int numCells = ms->numCells();
        for (int m = 0; m != numCells; ++m)
        {
            file.writeRow(m, units->otemperature(ms->gasTemperature(m)));
        }
    }
}

////////////////////////////////////////////////////////////////////
