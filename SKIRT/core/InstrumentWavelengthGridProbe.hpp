/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef WAVELENGTHGRIDPROBE_HPP
#define WAVELENGTHGRIDPROBE_HPP

#include "SpecialtyProbe.hpp"
class WavelengthGrid;

////////////////////////////////////////////////////////////////////

/** WavelengthGridProbe outputs a text column file with information on the wavelength grid for each
    instrument. Each file is named <tt>prefix_instr_wavelengths.txt</tt> so that it sits next to
    the files written by the corresponding instrument. For each wavelength bin, the file lists the
    characteristic wavelength, the wavelength bin width, and the left and right borders of the bin.
    */
class InstrumentWavelengthGridProbe : public SpecialtyProbe
{
    ITEM_CONCRETE(InstrumentWavelengthGridProbe, SpecialtyProbe, "wavelength grid: instruments")
        ATTRIBUTE_TYPE_DISPLAYED_IF(InstrumentWavelengthGridProbe, "Level2&Instrument")
    ITEM_END()

    //======================== Other Functions =======================

protected:
    /** This function performs probing. */
    void probe() override;

public:
    /** This function outputs a column text file for the specified wavelength grid in the format as
        described in the header of this class. It can be used from other probes to output
        wavelength grid details. */
    static void writeWavelengthGrid(Probe* item, const WavelengthGrid* wavelengthGrid, string filename,
                                    string description);
};

////////////////////////////////////////////////////////////////////

#endif
