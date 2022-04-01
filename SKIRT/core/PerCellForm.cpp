/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "PerCellForm.hpp"
#include "ProbeFormBridge.hpp"
#include "SpatialGrid.hpp"
#include "StringUtils.hpp"
#include "TextOutFile.hpp"

////////////////////////////////////////////////////////////////////

void PerCellForm::writeQuantity(const ProbeFormBridge* bridge) const
{
    // create a text column file and add the column definitions
    TextOutFile outfile(bridge->probe(), bridge->prefix(), bridge->description() + " per spatial cell");
    outfile.writeLine("# " + StringUtils::toUpperFirst(bridge->description()) + " per spatial cell");
    outfile.addColumn("spatial cell index", "", 'd');
    bridge->addColumnDefinitions(outfile);

    // allocate room for a single row
    int numValues = bridge->numValues();
    Array values(numValues);
    vector<double> columns(numValues + 1);

    // determine the number of spatial cells
    int numCells = bridge->grid() ? bridge->grid()->numCells() : 0;

    // write a line for each cell
    for (int m = 0; m != numCells; ++m)
    {
        bridge->valuesInCell(m, values);
        columns[0] = m;
        for (int p = 0; p != numValues; ++p) columns[p + 1] = values[p];
        outfile.writeRow(columns);
    }
}

////////////////////////////////////////////////////////////////////
