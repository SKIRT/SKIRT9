/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef TODDLERSSED_HPP
#define TODDLERSSED_HPP

#include "FamilySED.hpp"

////////////////////////////////////////////////////////////////////

/** A ToddlersSED class instance represents a single star-forming region SED taken from the
    TODDLERS template library, parameterized on age, metallicity, star formation efficiency, and
    cloud number density. See the ToddlersSEDFamily class for more information. */
class ToddlersSED : public FamilySED
{
    /** The enumeration type indicating the maximum PAH-to-dust fraction. */
    ENUM_DEF(PAHFraction, High, Low)
        ENUM_VAL(PAHFraction, High, "High PAH-to-dust fraction (4.6%)")
        ENUM_VAL(PAHFraction, Low, "Low PAH-to-dust fraction (1%)")
    ENUM_END()

    /** The enumeration type indicating the wavelength resolution. */
    ENUM_DEF(Resolution, Low, High)
        ENUM_VAL(Resolution, Low, "Low wavelength resolution (continuum and lines at R=300)")
        ENUM_VAL(Resolution, High, "High wavelength resolution (continuum at R=300 and lines at R=5e4)")
    ENUM_END()

    ITEM_CONCRETE(ToddlersSED, FamilySED, "a Toddlers SED for emission from star-forming regions")

        PROPERTY_ENUM(pahfraction, PAHFraction, "the maximum PAH-to-dust fraction")
        ATTRIBUTE_DEFAULT_VALUE(pahfraction, "High")

        PROPERTY_ENUM(resolution, Resolution, "the wavelength resolution")
        ATTRIBUTE_DEFAULT_VALUE(resolution, "Low")

        PROPERTY_DOUBLE(age, "system age")
        ATTRIBUTE_QUANTITY(age, "time")
        ATTRIBUTE_MIN_VALUE(age, "[0.1 Myr")
        ATTRIBUTE_MAX_VALUE(age, "30 Myr]")
        ATTRIBUTE_DEFAULT_VALUE(age, "2.5 Myr")

        PROPERTY_DOUBLE(metallicity, "system metallicity")
        ATTRIBUTE_MIN_VALUE(metallicity, "[0.001")
        ATTRIBUTE_MAX_VALUE(metallicity, "0.04]")
        ATTRIBUTE_DEFAULT_VALUE(metallicity, "0.02")

        PROPERTY_DOUBLE(SFE, "star formation efficiency ")
        ATTRIBUTE_MIN_VALUE(SFE, "[.01")
        ATTRIBUTE_MAX_VALUE(SFE, ".15]")
        ATTRIBUTE_DEFAULT_VALUE(SFE, ".025")

        PROPERTY_DOUBLE(cloudNumDensity, "natal cloud number density")
        ATTRIBUTE_QUANTITY(cloudNumDensity, "numbervolumedensity")
        ATTRIBUTE_MIN_VALUE(cloudNumDensity, "[10 /cm3")
        ATTRIBUTE_MAX_VALUE(cloudNumDensity, "2560 /cm3]")
        ATTRIBUTE_DEFAULT_VALUE(cloudNumDensity, "320 /cm3")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function returns a newly created SEDFamily object (which is already hooked into the
        simulation item hierachy so it will be automatically deleted) and stores the parameters for
        the specific %SED configured by the user in the specified array. */
    const SEDFamily* getFamilyAndParameters(Array& parameters) override;
};

////////////////////////////////////////////////////////////////////

#endif
