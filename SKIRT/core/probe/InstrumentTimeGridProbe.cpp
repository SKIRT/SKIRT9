/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "InstrumentTimeGridProbe.hpp"
#include "InstrumentSystem.hpp"
#include "TextOutFile.hpp"
#include "TimeGrid.hpp"
#include "TimeInstrument.hpp"
#include "Units.hpp"

////////////////////////////////////////////////////////////////////

void InstrumentTimeGridProbe::probe()
{
    // loop over instruments
    for (auto instrument : find<InstrumentSystem>()->instruments())
    {
        auto timeInstrument = dynamic_cast<TimeInstrument*>(instrument);
        if (timeInstrument)
        {
            auto units = find<Units>();

            // create a text file and add the columns
            string filename = timeInstrument->instrumentName() + "_timegrid";
            string description = "time grid for instrument " + timeInstrument->instrumentName();
            TextOutFile file(this, filename, description);
            file.addColumn("characteristic time", units->utimelag());
            file.addColumn("left border of time bin", units->utimelag());
            file.addColumn("right border of time bin", units->utimelag());
            file.addColumn("width of time bin", units->utimelag());

            // write the rows
            auto timeGrid = timeInstrument->timeGrid();
            int numBins = timeGrid->numBins();
            for (int k = 0; k != numBins; ++k)
            {
                double time = units->otimelag(timeGrid->time(k));
                double left = units->otimelag(timeGrid->left(k));
                double right = units->otimelag(timeGrid->right(k));
                double width = units->otimelag(timeGrid->width(k));
                file.writeRow(time, left, right, width);
            }
        }
    }
}

////////////////////////////////////////////////////////////////////
