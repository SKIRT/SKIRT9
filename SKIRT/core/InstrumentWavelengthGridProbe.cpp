/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "InstrumentWavelengthGridProbe.hpp"
#include "DisjointWavelengthGrid.hpp"
#include "Indices.hpp"
#include "Instrument.hpp"
#include "InstrumentSystem.hpp"
#include "TextOutFile.hpp"
#include "Units.hpp"
#include "WavelengthGrid.hpp"

////////////////////////////////////////////////////////////////////

void InstrumentWavelengthGridProbe::probe()
{
    // loop over instruments
    for (auto instrument : find<InstrumentSystem>()->instruments())
    {
        writeWavelengthGrid(this, instrument->instrumentWavelengthGrid(), instrument->instrumentName() + "_wavelengths",
                            "wavelengths for instrument " + instrument->instrumentName());
    }
}

////////////////////////////////////////////////////////////////////

void InstrumentWavelengthGridProbe::writeWavelengthGrid(Probe* item, const WavelengthGrid* wavelengthGrid,
                                                        string filename, string description)
{
    auto units = item->find<Units>();

    // create a text file and add the columns
    TextOutFile file(item, filename, description);
    file.addColumn("characteristic wavelength; " + units->swavelength(), units->uwavelength());
    file.addColumn("effective wavelength bin width; " + units->swavelength(), units->uwavelength());
    file.addColumn("left border of wavelength bin; " + units->swavelength(), units->uwavelength());
    file.addColumn("right border of wavelength bin; " + units->swavelength(), units->uwavelength());

    // write the rows
    for (int ell : Indices(wavelengthGrid->numBins(), units->rwavelength()))
    {
        double chara = units->owavelength(wavelengthGrid->wavelength(ell));
        double width = units->owavelength(wavelengthGrid->effectiveWidth(ell));
        double left = units->owavelength(wavelengthGrid->leftBorder(ell));
        double right = units->owavelength(wavelengthGrid->rightBorder(ell));

        // special treatment for output units with reverse ordering
        if (units->rwavelength())
        {
            // swap left and right borders
            std::swap(left, right);

            // fix the bin width
            if (dynamic_cast<const DisjointWavelengthGrid*>(wavelengthGrid))
            {
                // for wavelength grids with constant transmission across the bin,
                // the effective bin width is easily obtained from the borders
                width = right - left;
            }
            else
            {
                // for wavelength grids with variable transmission across the bin,
                // we apply the approximate correction factor dlambda^2 / lambda^2
                double lambda = wavelengthGrid->wavelength(ell);
                double dlambda = wavelengthGrid->effectiveWidth(ell);
                width *= (dlambda * dlambda) / (lambda * lambda);
            }
        }
        file.writeRow(chara, width, left, right);
    }
}

////////////////////////////////////////////////////////////////////
