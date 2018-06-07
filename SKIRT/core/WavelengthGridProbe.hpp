/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef WAVELENGTHGRIDPROBE_HPP
#define WAVELENGTHGRIDPROBE_HPP

#include "Probe.hpp"

////////////////////////////////////////////////////////////////////

/** WavelengthGridProbe outputs a text column file with information on the wavelength grid for each
    instrument. Each file is named <tt>prefix_instr_wavelengths.txt</tt> so that it sits next to
    the files written by the corresponding instrument. For each wavelength bin, the file lists the
    characteristic wavelength, the wavelength bin width, and the left and right borders of the bin.
    */
class WavelengthGridProbe : public Probe
{
    ITEM_CONCRETE(WavelengthGridProbe, Probe, "a probe of the instrument wavelength grids")
    ITEM_END()

    //======================== Other Functions =======================

public:
    /** This function performs probing after setup. */
    void probeSetup() override;
};

////////////////////////////////////////////////////////////////////

#endif
