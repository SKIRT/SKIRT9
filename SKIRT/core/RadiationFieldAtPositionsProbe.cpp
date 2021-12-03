/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "RadiationFieldAtPositionsProbe.hpp"
#include "Configuration.hpp"
#include "DisjointWavelengthGrid.hpp"
#include "Indices.hpp"
#include "InstrumentWavelengthGridProbe.hpp"
#include "MediumSystem.hpp"
#include "StringUtils.hpp"
#include "TextInFile.hpp"
#include "TextOutFile.hpp"
#include "Units.hpp"

////////////////////////////////////////////////////////////////////

void RadiationFieldAtPositionsProbe::probeRun()
{
    if (find<Configuration>()->hasRadiationField())
    {
        // output the mean intensity for each position in the input file
        {
            auto wavelengthGrid = find<Configuration>()->radiationFieldWLG();
            auto ms = find<MediumSystem>();
            auto grid = ms->grid();
            auto units = find<Units>();

            // open the input file
            TextInFile infile(this, filename(), "positions");
            infile.useColumns(useColumns());
            infile.addColumn("position x", "length", "pc");
            infile.addColumn("position y", "length", "pc");
            infile.addColumn("position z", "length", "pc");

            // create the output file
            TextOutFile outfile(this, itemName() + "_J", "mean intensity at positions");

            // write the header
            outfile.writeLine("# Mean radiation field intensities at imported positions");
            outfile.addColumn("position x", units->ulength());
            outfile.addColumn("position y", units->ulength());
            outfile.addColumn("position z", units->ulength());
            for (int ell : Indices(wavelengthGrid->numBins(), units->rwavelength()))
                outfile.addColumn(units->smeanintensity() + " at " + units->swavelength() + " = "
                                      + StringUtils::toString(units->owavelength(wavelengthGrid->wavelength(ell)), 'g')
                                      + " " + units->uwavelength(),
                                  units->umeanintensity());

            // write a line for each line in the input file
            double x, y, z;
            while (infile.readRow(x, y, z))
            {
                vector<double> values({units->olength(x), units->olength(y), units->olength(z)});
                int m = grid->cellIndex(Position(x, y, z));
                if (m >= 0)
                {
                    const Array& Jv = ms->meanIntensity(m);
                    for (int ell : Indices(wavelengthGrid->numBins(), units->rwavelength()))
                        values.push_back(units->omeanintensity(wavelengthGrid->wavelength(ell), Jv[ell]));
                }
                else
                {
                    for (int ell = 0; ell != wavelengthGrid->numBins(); ++ell) values.push_back(0.);
                }
                outfile.writeRow(values);
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
