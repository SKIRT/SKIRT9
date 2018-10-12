/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "DustAbsorptionPerCellProbe.hpp"
#include "Configuration.hpp"
#include "MediumSystem.hpp"
#include "SpatialGrid.hpp"
#include "StringUtils.hpp"
#include "TextOutFile.hpp"
#include "Units.hpp"
#include "WavelengthGrid.hpp"
#include "WavelengthGridProbe.hpp"

////////////////////////////////////////////////////////////////////

void DustAbsorptionPerCellProbe::probeRun()
{
    if (find<Configuration>()->hasRadiationField() && find<MediumSystem>()->hasDust())
    {
        // output the absorbed luminosity for each cell
        {
            auto wavelengthGrid = find<Configuration>()->radiationFieldWavelengthGrid();
            auto ms = find<MediumSystem>();
            auto grid = ms->grid();
            auto units = find<Units>();

            // create a text file
            TextOutFile file(this, itemName() + "_Labs", "dust absorption per cell");

            // write the header
            file.writeLine("# Spectral luminosity absorbed by dust per spatial cell");
            file.addColumn("spatial cell index", "", 'd');
            for (int ell=0; ell!=wavelengthGrid->numBins(); ++ell)
                file.addColumn(units->smonluminosity() + "^abs at lambda = "
                               + StringUtils::toString(units->owavelength(wavelengthGrid->wavelength(ell)), 'g')
                               + " " + units->uwavelength(), units->umonluminosity());

            // write a line for each cell
            int numCells = grid->numCells();
            for (int m=0; m!=numCells; ++m)
            {
                vector<double> values({ static_cast<double>(m) });
                const Array& Jv = ms->meanIntensity(m);
                double factor = 4.*M_PI * ms->volume(m);
                for (int ell=0; ell!=wavelengthGrid->numBins(); ++ell)
                {
                    double lambda = wavelengthGrid->wavelength(ell);
                    double Labs = Jv[ell] * factor * ms->opacityAbs(lambda, m, MaterialMix::MaterialType::Dust);
                    values.push_back(units->omonluminosityWavelength(lambda, Labs));
                }
                file.writeRow(values);
            }
        }

        // if requested, also output the wavelength grid
        if (writeWavelengthGrid())
        {
            WavelengthGridProbe::writeWavelengthGrid(this, find<Configuration>()->radiationFieldWavelengthGrid(),
                                                     itemName() + "_wavelengths", "wavelengths for mean intensity");
        }
    }
}

////////////////////////////////////////////////////////////////////
