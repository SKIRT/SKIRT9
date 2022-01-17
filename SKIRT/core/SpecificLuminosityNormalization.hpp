/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef SPECIFICLUMINOSITYNORMALIZATION_HPP
#define SPECIFICLUMINOSITYNORMALIZATION_HPP

#include "LuminosityNormalization.hpp"

////////////////////////////////////////////////////////////////////

/** A SpecificLuminosityNormalization instance sets the normalization of a primary source by
    specifying the specific luminosity (radiative power per units of wavelength, frequency, or
    energy) at a certain wavelength. This normalization cannot be used with LineSED subclasses, for
    which the specific luminosity is undefined. An attempt to do so will result in a fatal error.
    */
class SpecificLuminosityNormalization : public LuminosityNormalization
{
    /** The enumeration type indicating the specific luminosity unit style, e.g. whether to use
        specific luminosity per unit of wavelength, frequency or energy. */
    ENUM_DEF(UnitStyle, neutralmonluminosity, wavelengthmonluminosity, frequencymonluminosity, energymonluminosity)
        ENUM_VAL(UnitStyle, neutralmonluminosity, "neutral: λ L_λ = ν L_ν")
        ENUM_VAL(UnitStyle, wavelengthmonluminosity, "per unit of wavelength: L_λ")
        ENUM_VAL(UnitStyle, frequencymonluminosity, "per unit of frequency: L_ν")
        ENUM_VAL(UnitStyle, energymonluminosity, "counts per unit of energy: L_E")
    ENUM_END()

    ITEM_CONCRETE(SpecificLuminosityNormalization, LuminosityNormalization,
                  "source normalization through the specific luminosity at a given wavelength")
        ATTRIBUTE_TYPE_ALLOWED_IF(SpecificLuminosityNormalization, "ContSED")

        PROPERTY_DOUBLE(wavelength, "the wavelength at which to provide the specific luminosity")
        ATTRIBUTE_QUANTITY(wavelength, "wavelength")
        ATTRIBUTE_MIN_VALUE(wavelength, "1 pm")
        ATTRIBUTE_MAX_VALUE(wavelength, "1 m")

        PROPERTY_ENUM(unitStyle, UnitStyle, "the luminosity unit style")
        ATTRIBUTE_DEFAULT_VALUE(unitStyle, "wavelengthmonluminosity")

        PROPERTY_DOUBLE(specificLuminosity, "the specific luminosity at the given wavelength")
        ATTRIBUTE_QUANTITY(specificLuminosity, "@unitStyle")
        ATTRIBUTE_MIN_VALUE(specificLuminosity, "]0")

    ITEM_END()

    //======================== Other Functions =======================

public:
    /** This function returns the luminosity of a source with the specified %SED, limited to the
        source's wavelength range. Given that the luminosity of the spectrum described by the SED
        object over the source wavelength range is normalized to unity, the requested luminosity is
        obtained by dividing the user-configured specific luminosity by the normalized specific
        luminosity in the %SED at the same wavelength. */
    double luminosityForSED(SED* sed) const override;
};

////////////////////////////////////////////////////////////////////

#endif
