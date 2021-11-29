/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "ElectronTemperaturePerCellProbe.hpp"
#include "Configuration.hpp"
#include "MediumSystem.hpp"
#include "TextOutFile.hpp"
#include "Units.hpp"

////////////////////////////////////////////////////////////////////

void ElectronTemperaturePerCellProbe::probeSetup()
{
    if (find<Configuration>()->hasMedium() && find<MediumSystem>()->hasElectrons())
    {
        // locate the medium system
        auto ms = find<MediumSystem>();
        auto units = find<Units>();

        // create a text file
        TextOutFile file(this, itemName() + "_T", "electron temperature per cell");

        // write the header
        file.writeLine("# Electron temperature per spatial cell");
        file.addColumn("spatial cell index", "", 'd');
        file.addColumn("electron temperature", units->utemperature(), 'g');

        // write a line for each cell
        int numCells = ms->numCells();
        for (int m = 0; m != numCells; ++m)
        {
            file.writeRow(m, units->otemperature(ms->indicativeElectronTemperature(m)));
        }
    }
}

////////////////////////////////////////////////////////////////////
