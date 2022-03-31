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
    wavelength-dependent material properties may be required for the configured wavelength.

    In addition, this class offers an option for the user to decide whether the probe should be
    performed after setup or after the full simulation run. This functionality is unrelated to the
    wavelength property but (today) it happens to be required by the same subclasses so it is
    implemented here for convenience. */
class AbstractWavelengthProbe : public Probe, public MaterialWavelengthRangeInterface
{
    /** The enumeration type indicating when probing occurs. */
    ENUM_DEF(ProbeAfter, Setup, Run)
        ENUM_VAL(ProbeAfter, Setup, "after setup")
        ENUM_VAL(ProbeAfter, Run, "after the complete simulation run")
    ENUM_END()

    ITEM_ABSTRACT(AbstractWavelengthProbe, Probe, "a probe requiring a wavelength value")

        PROPERTY_DOUBLE(wavelength, "the wavelength at which to determine the optical depth")
        ATTRIBUTE_QUANTITY(wavelength, "wavelength")
        ATTRIBUTE_MIN_VALUE(wavelength, "1 pm")
        ATTRIBUTE_MAX_VALUE(wavelength, "1 m")
        ATTRIBUTE_DEFAULT_VALUE(wavelength, "0.55 micron")
        ATTRIBUTE_DISPLAYED_IF(wavelength, "Level2")

        ATTRIBUTE_SUB_PROPERTIES_HERE()

        PROPERTY_ENUM(probeAfter, ProbeAfter, "perform the probe after")
        ATTRIBUTE_DEFAULT_VALUE(probeAfter, "Setup")
        ATTRIBUTE_DISPLAYED_IF(probeAfter, "DynamicState")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function returns an enumeration indicating when probing for this probe should be
        performed corresponding to the configured value of the \em probeAfter property. */
    When when() const override;

    //======================== Other Functions =======================

public:
    /** This function returns a wavelength range corresponding to the user-configured wavelength,
        indicating that wavelength-dependent material properties may be required for this
        wavelength. */
    Range wavelengthRange() const override;
};

////////////////////////////////////////////////////////////////////

#endif
