/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef LUMINOSITYPROBE_HPP
#define LUMINOSITYPROBE_HPP

#include "Probe.hpp"
#include "WavelengthGrid.hpp"

////////////////////////////////////////////////////////////////////

/** LuminosityProbe outputs a text column file with the primary source luminosities on a specified
    wavelength grid (or on the default instrument wavelength grid). The file is named
    <tt>prefix_luminosities.txt</tt>. The first three columns list the characteristic wavelength of
    the wavelength bin, the total specific luminosity at that wavelength, and the total luminosity
    integrated over the bin (in both cases summed over all primary sources). Furthermore, there is
    an additional column for each primary source, listing the fractional contribution of that
    source. The columns are in the same order as the sources appear in the configuration file. If
    there is only one source, there is a single column containing one's.

    Note that, by definition, the luminosity is zero outside of the primary source wavelength range
    configured in the source system. */
class LuminosityProbe : public Probe
{
    ITEM_CONCRETE(LuminosityProbe, Probe, "the primary source luminosities")

    PROPERTY_ITEM(wavelengthGrid, WavelengthGrid, "the wavelength grid for the luminosity probe")
        ATTRIBUTE_REQUIRED_IF(wavelengthGrid, "!DefaultInstrumentWavelengthGrid")

    ITEM_END()

    //======================== Other Functions =======================

public:
    /** This function performs probing after setup. */
    void probeSetup() override;
};

////////////////////////////////////////////////////////////////////

#endif
