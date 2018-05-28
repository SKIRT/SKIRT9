/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef WAVELENGTHPROBE_HPP
#define WAVELENGTHPROBE_HPP

#include "Probe.hpp"

////////////////////////////////////////////////////////////////////

/** WavelengthProbe outputs a text column file with the wavelength grid for each instrument. */
class WavelengthProbe : public Probe
{
    ITEM_CONCRETE(WavelengthProbe, Probe, "a probe of the instrument wavelength grids")
    ITEM_END()

    //======================== Other Functions =======================

public:
    /** This function performs probing after setup. */
    void probeSetup() override;
};

////////////////////////////////////////////////////////////////////

#endif
