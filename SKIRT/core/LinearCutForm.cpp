/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "LinearCutForm.hpp"
#include "ProbeFormBridge.hpp"
#include "StringUtils.hpp"
#include "TextOutFile.hpp"
#include "Units.hpp"

////////////////////////////////////////////////////////////////////

void LinearCutForm::writeQuantity(const ProbeFormBridge* bridge) const
{
    auto units = bridge->units();

    // create a text column file and add the column definitions
    TextOutFile outfile(bridge->probe(), bridge->prefix(), bridge->description() + " along a line segment");
    outfile.writeLine("# " + StringUtils::toUpperFirst(bridge->description()) + " along a line segment");
    outfile.addColumn("distance from starting point", units->ulength());
    bridge->addColumnDefinitions(outfile);

    // allocate room for a single row
    int numValues = bridge->numValues();
    Array values(numValues);
    vector<double> columns(numValues + 1);

    // get the line segment start and end points
    Position p1(_startX, _startY, _startZ);
    Position p2(_endX, _endY, _endZ);

    // write a line for each sample
    for (int i = 0; i != _numSamples; ++i)
    {
        // determine the sample position and the distance along the line
        double fraction = static_cast<double>(i) / static_cast<double>(_numSamples - 1);
        Position p(p1 + fraction * (p2 - p1));
        double distance = (p - p1).norm();

        // get the corresponding quantity value and write the row
        bridge->valuesAtPosition(p, values);
        columns[0] = units->olength(distance);
        for (int p = 0; p != numValues; ++p) columns[p + 1] = values[p];
        outfile.writeRow(columns);
    }
}

////////////////////////////////////////////////////////////////////
