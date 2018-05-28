/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "WavelengthProbe.hpp"
#include "Instrument.hpp"
#include "InstrumentSystem.hpp"
#include "TextOutFile.hpp"
#include "Units.hpp"
#include "WavelengthGrid.hpp"

////////////////////////////////////////////////////////////////////

void WavelengthProbe::probeSetup()
{
    auto units = find<Units>();

    // loop over instruments
    for (auto instrument : find<InstrumentSystem>()->instruments())
    {
        auto wavelengthGrid = instrument->instrumentWavelengthGrid();
        auto instrumentName = instrument->instrumentName();

        // create a text file and add the columns
        TextOutFile file(this, instrumentName+"_wavelengths", "wavelengths for instrument " + instrumentName);
        file.addColumn("characteristic wavelength (" + units->uwavelength() + ")", 'e', 9);
        file.addColumn("wavelength bin width (" + units->uwavelength() + ")", 'e', 9);
        file.addColumn("left border of wavelength bin (" + units->uwavelength() + ")", 'e', 9);
        file.addColumn("right border of wavelength bin (" + units->uwavelength() + ")", 'e', 9);

        // write the rows
        int numWavelengths = wavelengthGrid->numWavelengths();
        for (int ell=0; ell!=numWavelengths; ++ell)
        {
            file.writeRow(units->owavelength(wavelengthGrid->lambda(ell)),
                          units->owavelength(wavelengthGrid->dlambda(ell)),
                          units->owavelength(wavelengthGrid->lambdaLeft(ell)),
                          units->owavelength(wavelengthGrid->lambdaRight(ell)) );
        }
    }
}

////////////////////////////////////////////////////////////////////
