/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "WavelengthGridProbe.hpp"
#include "Instrument.hpp"
#include "InstrumentSystem.hpp"
#include "TextOutFile.hpp"
#include "Units.hpp"
#include "WavelengthGrid.hpp"

////////////////////////////////////////////////////////////////////

void WavelengthGridProbe::probeSetup()
{
    auto units = find<Units>();

    // loop over instruments
    for (auto instrument : find<InstrumentSystem>()->instruments())
    {
        auto wavelengthGrid = instrument->instrumentWavelengthGrid();
        auto instrumentName = instrument->instrumentName();

        // create a text file and add the columns
        TextOutFile file(this, instrumentName+"_wavelengths", "wavelengths for instrument " + instrumentName);
        file.addColumn("characteristic wavelength", units->uwavelength());
        file.addColumn("effective wavelength bin width", units->uwavelength());
        file.addColumn("left border of wavelength bin", units->uwavelength());
        file.addColumn("right border of wavelength bin", units->uwavelength());

        // write the rows
        int numWavelengths = wavelengthGrid->numBins();
        for (int ell=0; ell!=numWavelengths; ++ell)
        {
            file.writeRow(units->owavelength(wavelengthGrid->wavelength(ell)),
                          units->owavelength(wavelengthGrid->effectiveWidth(ell)),
                          units->owavelength(wavelengthGrid->leftBorder(ell)),
                          units->owavelength(wavelengthGrid->rightBorder(ell)) );
        }
    }
}

////////////////////////////////////////////////////////////////////
