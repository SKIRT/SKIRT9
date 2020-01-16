/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "MediumVelocityPerCellProbe.hpp"
#include "Configuration.hpp"
#include "MediumSystem.hpp"
#include "SpatialGrid.hpp"
#include "StringUtils.hpp"
#include "TextOutFile.hpp"
#include "Units.hpp"
#include "WavelengthGrid.hpp"

////////////////////////////////////////////////////////////////////

void MediumVelocityPerCellProbe::probeSetup()
{
    if (find<Configuration>()->hasMovingMedia())
    {
        // locate the medium system
        auto ms = find<MediumSystem>();
        auto units = find<Units>();

        // create a text file
        TextOutFile file(this, itemName() + "_v", "medium velocity per cell");

        // write the header
        file.writeLine("# Medium velocity per spatial cell");
        file.addColumn("spatial cell index", "", 'd');
        file.addColumn("x component of velocity", units->uvelocity());
        file.addColumn("y component of velocity", units->uvelocity());
        file.addColumn("z component of velocity", units->uvelocity());

        // write a line for each cell
        int numCells = ms->numCells();
        for (int m = 0; m != numCells; ++m)
        {
            Vec v = ms->bulkVelocity(m);
            file.writeRow(m, units->omagneticfield(v.x()), units->omagneticfield(v.y()), units->omagneticfield(v.z()));
        }
    }
}

////////////////////////////////////////////////////////////////////
