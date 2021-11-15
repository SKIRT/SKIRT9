/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef ABSTRACTWAVELENGTHPROBE_HPP
#define ABSTRACTWAVELENGTHPROBE_HPP

#include "MaterialWavelengthRangeInterface.hpp"
#include "Probe.hpp"

////////////////////////////////////////////////////////////////////

/** AbstractWavelengthProbe is a base class for probes that require a configurable wavelength
    property. It implements the MaterialWavelengthRangeInterface to indicate that
    wavelength-dependent material properties may be required for the configured wavelength. */
class AbstractWavelengthProbe : public Probe, public MaterialWavelengthRangeInterface
{
    ITEM_ABSTRACT(AbstractWavelengthProbe, Probe, "a probe requiring a wavelength value")

        PROPERTY_DOUBLE(wavelength, "the wavelength at which to determine the optical depth")
        ATTRIBUTE_QUANTITY(wavelength, "wavelength")
        ATTRIBUTE_MIN_VALUE(wavelength, "1 pm")
        ATTRIBUTE_MAX_VALUE(wavelength, "1 m")
        ATTRIBUTE_DEFAULT_VALUE(wavelength, "0.55 micron")
        ATTRIBUTE_DISPLAYED_IF(wavelength, "Level2")

    ITEM_END()

    //======================== Other Functions =======================

public:
    /** This function returns a wavelength range corresponding to the user-configured wavelength,
        indicating that wavelength-dependent material properties may be required for this
        wavelength. */
    Range wavelengthRange() const override;
};

////////////////////////////////////////////////////////////////////

#endif
