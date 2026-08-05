/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "AtPositionsForm.hpp"
#include "ProbeFormBridge.hpp"
#include "StringUtils.hpp"
#include "TextInFile.hpp"
#include "TextOutFile.hpp"
#include "Units.hpp"

////////////////////////////////////////////////////////////////////

void AtPositionsForm::writeQuantity(const ProbeFormBridge* bridge) const
{
    auto units = bridge->units();

    // create a text column file and add the column definitions
    TextOutFile outfile(bridge->probe(), bridge->prefix(), bridge->description() + " at imported positions");
    outfile.writeLine("# " + StringUtils::toUpperFirst(bridge->description()) + " at imported positions");
    outfile.addColumn("position x", units->ulength());
    outfile.addColumn("position y", units->ulength());
    outfile.addColumn("position z", units->ulength());
    bridge->addColumnDefinitions(outfile);

    // allocate room for a single row
    int numValues = bridge->numValues();
    Array values(numValues);
    vector<double> columns(numValues + 3);

    // open the input file
    TextInFile infile(this, filename(), "positions");
    infile.useColumns(useColumns());
    infile.addColumn("position x", "length", "pc");
    infile.addColumn("position y", "length", "pc");
    infile.addColumn("position z", "length", "pc");

    // write a line for each line in the input file
    double x, y, z;
    while (infile.readRow(x, y, z))
    {
        bridge->valuesAtPosition(Position(x, y, z), values);
        columns[0] = units->olength(x);
        columns[1] = units->olength(y);
        columns[2] = units->olength(z);
        for (int p = 0; p != numValues; ++p) columns[p + 3] = values[p];
        outfile.writeRow(columns);
    }
}

////////////////////////////////////////////////////////////////////
