/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "InstrumentWavelengthGridProbe.hpp"
#include "Indices.hpp"
#include "Instrument.hpp"
#include "InstrumentSystem.hpp"
#include "TextOutFile.hpp"
#include "Units.hpp"
#include "WavelengthGrid.hpp"

////////////////////////////////////////////////////////////////////

void InstrumentWavelengthGridProbe::probeSetup()
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
        file.writeRow(units->owavelength(wavelengthGrid->wavelength(ell)),
                      units->owavelength(wavelengthGrid->effectiveWidth(ell)),
                      units->owavelength(wavelengthGrid->leftBorder(ell)),
                      units->owavelength(wavelengthGrid->rightBorder(ell)));
    }
}

////////////////////////////////////////////////////////////////////
