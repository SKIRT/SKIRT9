/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "MeridionalDustTemperatureCutProbe.hpp"
#include "Configuration.hpp"
#include "FatalError.hpp"
#include "MediumSystem.hpp"
#include "StringUtils.hpp"
#include "TextOutFile.hpp"
#include "Units.hpp"
#include "WavelengthGrid.hpp"

////////////////////////////////////////////////////////////////////

void MeridionalDustTemperatureCutProbe::probeRun()
{
    if (find<Configuration>()->hasRadiationField() && find<MediumSystem>()->hasDust())
    {
        // locate the medium system
        auto ms = find<MediumSystem>();
        auto units = find<Units>();

        // create a text file
        TextOutFile file(this, itemName() + "_T", "dust temperature along a line segment");

        // write the header
        file.writeLine("# Indicative dust temperature along a meridian");
        file.addColumn("inclination", units->uposangle());
        file.addColumn("indicative dust temperature", units->utemperature(), 'g');

        // write a line for each sample
        for (int i=0; i!=_numSamples; ++i)
        {
            // determine the sample inclination and position
            double fraction = static_cast<double>(i) / static_cast<double>(_numSamples-1);
            double inclination = fraction * M_PI;
            Position p(_radius * Direction(inclination, _azimuth));

            // calculate the corresponding indicative dust temperature and write the row
            file.writeRow(units->oposangle(inclination), units->otemperature(ms->indicativeDustTemperature(p)));
        }
    }
}

////////////////////////////////////////////////////////////////////
