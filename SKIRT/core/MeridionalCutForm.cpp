/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "MeridionalCutForm.hpp"
#include "ProbeFormBridge.hpp"
#include "StringUtils.hpp"
#include "TextOutFile.hpp"
#include "Units.hpp"

////////////////////////////////////////////////////////////////////

void MeridionalCutForm::writeQuantity(const ProbeFormBridge* bridge) const
{
    auto units = bridge->units();

    // create a text column file and add the column definitions
    TextOutFile outfile(bridge->probe(), bridge->prefix(), bridge->description() + " along a meridian");
    outfile.writeLine("# " + StringUtils::toUpperFirst(bridge->description()) + " along a meridian");
    outfile.addColumn("inclination", units->uposangle());
    bridge->addColumnDefinitions(outfile);

    // allocate room for a single row
    int numValues = bridge->numValues();
    Array values(numValues);
    vector<double> columns(numValues + 1);

    // write a line for each sample
    for (int i = 0; i != _numSamples; ++i)
    {
        // determine the sample inclination and position
        double fraction = static_cast<double>(i) / static_cast<double>(_numSamples - 1);
        double inclination = fraction * M_PI;
        Position p(_radius * Direction(inclination, _azimuth));

        // get the corresponding quantity value and write the row
        bridge->valuesAtPosition(p, values);
        columns[0] = units->oposangle(inclination);
        for (int p = 0; p != numValues; ++p) columns[p + 1] = values[p];
        outfile.writeRow(columns);
    }
}

////////////////////////////////////////////////////////////////////
