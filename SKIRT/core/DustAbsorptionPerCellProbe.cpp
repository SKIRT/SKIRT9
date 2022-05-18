/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "DustAbsorptionPerCellProbe.hpp"
#include "Configuration.hpp"
#include "DisjointWavelengthGrid.hpp"
#include "Indices.hpp"
#include "InstrumentWavelengthGridProbe.hpp"
#include "MediumSystem.hpp"
#include "StringUtils.hpp"
#include "TextOutFile.hpp"
#include "Units.hpp"

////////////////////////////////////////////////////////////////////

Probe::When DustAbsorptionPerCellProbe::when() const
{
    return When::Run;
}

////////////////////////////////////////////////////////////////////

void DustAbsorptionPerCellProbe::probe()
{
    if (find<Configuration>()->hasRadiationField() && find<MediumSystem>()->hasDust())
    {
        // output the absorbed luminosity for each cell
        {
            auto wavelengthGrid = find<Configuration>()->radiationFieldWLG();
            auto ms = find<MediumSystem>();
            auto grid = ms->grid();
            auto units = find<Units>();

            // create a text file
            TextOutFile file(this, itemName() + "_Labs", "dust absorption per cell");

            // write the header
            file.writeLine("# Spectral luminosity absorbed by dust per spatial cell");
            file.addColumn("spatial cell index", "", 'd');
            for (int ell : Indices(wavelengthGrid->numBins(), units->rwavelength()))
                file.addColumn(units->smonluminosity() + "^abs at " + units->swavelength() + " = "
                                   + StringUtils::toString(units->owavelength(wavelengthGrid->wavelength(ell)), 'g')
                                   + " " + units->uwavelength(),
                               units->umonluminosity());

            // write a line for each cell
            int numCells = grid->numCells();
            for (int m = 0; m != numCells; ++m)
            {
                vector<double> values({static_cast<double>(m)});
                const Array& Jv = ms->meanIntensity(m);
                double factor = 4. * M_PI * ms->volume(m);
                for (int ell : Indices(wavelengthGrid->numBins(), units->rwavelength()))
                {
                    double lambda = wavelengthGrid->wavelength(ell);
                    double Labs = Jv[ell] * factor * ms->opacityAbs(lambda, m, MaterialMix::MaterialType::Dust);
                    values.push_back(units->omonluminosity(lambda, Labs));
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
