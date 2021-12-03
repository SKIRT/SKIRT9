/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "RadiationFieldPerCellProbe.hpp"
#include "Configuration.hpp"
#include "DisjointWavelengthGrid.hpp"
#include "InstrumentWavelengthGridProbe.hpp"
#include "MediumSystem.hpp"
#include "StringUtils.hpp"
#include "TextOutFile.hpp"
#include "Units.hpp"

////////////////////////////////////////////////////////////////////

void RadiationFieldPerCellProbe::probeRun()
{
    if (find<Configuration>()->hasRadiationField())
    {
        // output the mean intensity for each cell
        {
            auto wavelengthGrid = find<Configuration>()->radiationFieldWLG();
            auto ms = find<MediumSystem>();
            auto grid = ms->grid();
            auto units = find<Units>();

            // create a text file
            TextOutFile file(this, itemName() + "_J", "mean intensity per cell");

            // write the header
            file.writeLine("# Mean radiation field intensities per spatial cell");
            file.addColumn("spatial cell index", "", 'd');
            for (int ell = 0; ell != wavelengthGrid->numBins(); ++ell)
                file.addColumn(units->smeanintensity() + " at " + units->swavelength() + " = "
                                   + StringUtils::toString(units->owavelength(wavelengthGrid->wavelength(ell)), 'g')
                                   + " " + units->uwavelength(),
                               units->umeanintensity());

            // write a line for each cell
            int numCells = grid->numCells();
            for (int m = 0; m != numCells; ++m)
            {
                vector<double> values({static_cast<double>(m)});
                const Array& Jv = ms->meanIntensity(m);
                for (int ell = 0; ell != wavelengthGrid->numBins(); ++ell)
                {
                    values.push_back(units->omeanintensity(wavelengthGrid->wavelength(ell), Jv[ell]));
                }
                file.writeRow(values);
            }
        }

        // if requested, also output the wavelength grid
        if (writeWavelengthGrid())
        {
            InstrumentWavelengthGridProbe::writeWavelengthGrid(this, find<Configuration>()->radiationFieldWLG(),
                                                               itemName() + "_wavelengths",
                                                               "wavelengths for mean intensity");
        }
    }
}

////////////////////////////////////////////////////////////////////
