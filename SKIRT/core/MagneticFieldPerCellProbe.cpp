/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "MagneticFieldPerCellProbe.hpp"
#include "Configuration.hpp"
#include "MediumSystem.hpp"
#include "SpatialGrid.hpp"
#include "StringUtils.hpp"
#include "TextOutFile.hpp"
#include "Units.hpp"
#include "WavelengthGrid.hpp"

////////////////////////////////////////////////////////////////////

void MagneticFieldPerCellProbe::probeRun()
{
    if (find<Configuration>()->hasPanRadiationField() && find<MediumSystem>()->hasDust())
    {
        // locate the medium system
        auto ms = find<MediumSystem>();
        auto units = find<Units>();

        // create a text file
        TextOutFile file(this, itemName() + "_B", "magnetic field per cell");

        // write the header
        file.writeLine("# Magnetic field per spatial cell");
        file.addColumn("spatial cell index", "", 'd');
        file.addColumn("x component of magnetic field", units->umagneticfield(), 'g');
        file.addColumn("y component of magnetic field", units->umagneticfield(), 'g');
        file.addColumn("z component of magnetic field", units->umagneticfield(), 'g');

        // write a line for each cell
        int numCells = ms->grid()->numCells();
        for (int m = 0; m != numCells; ++m)
        {
            Vec magneticField = ms->magneticField(m);
            file.writeRow(m, units->omagneticfield(magneticField.x()), units->omagneticfield(magneticField.y()),
                          units->omagneticfield(magneticField.z()));
        }
    }
}

////////////////////////////////////////////////////////////////////
