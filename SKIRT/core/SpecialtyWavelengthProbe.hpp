/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef SPECIALTYWAVELENGTHPROBE_HPP
#define SPECIALTYWAVELENGTHPROBE_HPP

#include "MaterialWavelengthRangeInterface.hpp"
#include "SpecialtyWhenProbe.hpp"

////////////////////////////////////////////////////////////////////

/** SpecialtyWavelengthProbe is a base class for probes that require a configurable wavelength
    property. It implements the MaterialWavelengthRangeInterface to indicate that
    wavelength-dependent material properties may be required for the configured wavelength.

    This class derives from the SpecialtyWhenProbe class because all subclasses happen to need the
    option for the user to decide whether the probe should be performed after setup or after the
    full simulation run. */
class SpecialtyWavelengthProbe : public SpecialtyWhenProbe, public MaterialWavelengthRangeInterface
{
    ITEM_ABSTRACT(SpecialtyWavelengthProbe, SpecialtyWhenProbe, "a specialty wavelength probe")

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
