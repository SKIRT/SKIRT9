/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef ABSTRACTWAVELENGTHGRIDPROBE_HPP
#define ABSTRACTWAVELENGTHGRIDPROBE_HPP

#include "Probe.hpp"
#include "WavelengthGrid.hpp"

////////////////////////////////////////////////////////////////////

/** AbstractWavelengthGridProbe is a base class for probes that require a configurable wavelength
    grid property. Deriving such probes from a common base class allows accessing the wavelength
    grid through the base class interface. */
class AbstractWavelengthGridProbe : public Probe
{
    ITEM_CONCRETE(AbstractWavelengthGridProbe, Probe, "a probe requiring a wavelength grid")

    PROPERTY_ITEM(wavelengthGrid, WavelengthGrid, "the wavelength grid for this probe")
        ATTRIBUTE_RELEVANT_IF(wavelengthGrid, "Panchromatic")
        ATTRIBUTE_REQUIRED_IF(wavelengthGrid, "!DefaultInstrumentWavelengthGrid")
        ATTRIBUTE_DISPLAYED_IF(wavelengthGrid, "Level2")

    ITEM_END()
};

////////////////////////////////////////////////////////////////////

#endif
