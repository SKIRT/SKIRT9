/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef LISTSED_HPP
#define LISTSED_HPP

#include "TabulatedSED.hpp"

////////////////////////////////////////////////////////////////////

/** A ListSED object represents a spectral energy distribution that is fully specified inside the
    configuration file (i.e. without referring to an input file). It is intended for use in cases
    where there are just a few wavelength/luminosity pairs, but nothing keeps the user from
    specifying a long list.

    The specific luminosity values can in fact be given in arbitrary units because the %SED will be
    normalized after being loaded. The luminosity outside the range indicated by the first and the
    last wavelength in the list is considered to be zero. */
class ListSED : public TabulatedSED
{
    /** The enumeration type indicating the specific luminosity unit style, e.g. whether to use
        specific luminosity per unit of wavelength or per unit of frequency. */
    ENUM_DEF(UnitStyle, wavelengthmonluminosity, frequencymonluminosity, neutralmonluminosity)
    ENUM_VAL(UnitStyle, wavelengthmonluminosity, "per unit of wavelength: L_λ")
    ENUM_VAL(UnitStyle, frequencymonluminosity, "per unit of frequency: L_ν")
    ENUM_VAL(UnitStyle, neutralmonluminosity, "neutral: λ L_λ = ν L_ν")
    ENUM_END()

    ITEM_CONCRETE(ListSED, TabulatedSED, "a spectral energy distribution specified inside the configuration file")

    PROPERTY_DOUBLE_LIST(wavelengths, "the wavelengths at which to specify the specific luminosity")
        ATTRIBUTE_QUANTITY(wavelengths, "wavelength")
        ATTRIBUTE_MIN_VALUE(wavelengths, "1 A")
        ATTRIBUTE_MAX_VALUE(wavelengths, "1 m")

    PROPERTY_ENUM(unitStyle, UnitStyle, "the luminosity unit style")
        ATTRIBUTE_DEFAULT_VALUE(unitStyle, "wavelengthmonluminosity")

    PROPERTY_DOUBLE_LIST(specificLuminosities, "the specific luminosities at each of the given wavelengths")
        ATTRIBUTE_QUANTITY(specificLuminosities, "@unitStyle")
        ATTRIBUTE_MIN_VALUE(specificLuminosities, "]0")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function verifies that the number of configured wavelengths and luminosities
        match, converts the luminosities to per-wavelength style if needed, and loads the data
        into the specified arrays. */
    void getWavelengthsAndLuminosities(Array& lambdav, Array& pv) const override;
};

////////////////////////////////////////////////////////////////////

#endif
