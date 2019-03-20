/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef ABSTRACTWAVELENGTHPROBE_HPP
#define ABSTRACTWAVELENGTHPROBE_HPP

#include "Probe.hpp"

////////////////////////////////////////////////////////////////////

/** AbstractWavelengthProbe is a base class for probes that require a configurable wavelength
    property. Deriving such probes from a common base class allows accessing the wavelength value
    through the base class interface. */
class AbstractWavelengthProbe : public Probe
{
    ITEM_CONCRETE(AbstractWavelengthProbe, Probe, "a probe requiring a wavelength value")
        ATTRIBUTE_TYPE_DISPLAYED_IF(SpatialGridConvergenceProbe, "Medium&SpatialGrid")

    PROPERTY_DOUBLE(wavelength, "the wavelength at which to determine the optical depth")
        ATTRIBUTE_QUANTITY(wavelength, "wavelength")
        ATTRIBUTE_MIN_VALUE(wavelength, "1 Angstrom")
        ATTRIBUTE_MAX_VALUE(wavelength, "1 m")
        ATTRIBUTE_DEFAULT_VALUE(wavelength, "0.55 micron")
        ATTRIBUTE_DISPLAYED_IF(wavelength, "Level2")

    ITEM_END()
};

////////////////////////////////////////////////////////////////////

#endif
