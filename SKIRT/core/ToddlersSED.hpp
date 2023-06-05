/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef TODDLERSSED_HPP
#define TODDLERSSED_HPP

#include "FamilySED.hpp"

////////////////////////////////////////////////////////////////////

/** ToddlersSED is a class representing the family of 
    Toddlers star-forming regions template SEDs, parameterized on age, metallicity,
    star formation efficiency, cloud number density and scaled by particle mass. 
    TODDLERS = Time evolution of Dust Diagnostics and Line Emission from Regions
    containing young Stars.
 */
class ToddlersSED : public FamilySED
{
    /** The enumeration type indicating the maximum PAH to Dust fraction value. */
    ENUM_DEF(PAHfraction, High, Low)
        ENUM_VAL(PAHfraction, High, "High PAH fraction")
        ENUM_VAL(PAHfraction, Low, "Low PAH fraction")
    ENUM_END()

    /** The enumeration type indicating the wavelength resolution. */
    ENUM_DEF(Resolution, Low, High)
        ENUM_VAL(Resolution, Low, "Low wavelength resolution (Cloudy Default, R=300)")
        ENUM_VAL(Resolution, High, "High wavelength resolution (Lines at R=1e5, Continua at R=300)")
    ENUM_END()

    ITEM_CONCRETE(ToddlersSED, FamilySED, "a Toddlers SED family for emission from star-forming regions")

        PROPERTY_ENUM(pahfraction, PAHfraction, "the maximum PAH to Dust fraction value")
        ATTRIBUTE_DEFAULT_VALUE(pahfraction, "High")
        PROPERTY_ENUM(resolution, Resolution,  "the wavelength resolution")
        ATTRIBUTE_DEFAULT_VALUE(resolution, "low")

        PROPERTY_DOUBLE(age, "Hii region age")
        ATTRIBUTE_QUANTITY(age, "time")
        ATTRIBUTE_MIN_VALUE(age, "[0.1 Myr")
        ATTRIBUTE_MAX_VALUE(age, "30 Myr]")
        ATTRIBUTE_DEFAULT_VALUE(age, "2.5 Myr")

        PROPERTY_DOUBLE(metallicity, "Hii region metallicity")
        ATTRIBUTE_MIN_VALUE(metallicity, "[0.001")
        ATTRIBUTE_MAX_VALUE(metallicity, "0.04]")
        ATTRIBUTE_DEFAULT_VALUE(metallicity, "0.02")

        PROPERTY_DOUBLE(SFE, "Initial star formation efficiency M*/Mtot ")
        ATTRIBUTE_MIN_VALUE(SFE, "[.01")
        ATTRIBUTE_MAX_VALUE(SFE, ".15]")
        ATTRIBUTE_DEFAULT_VALUE(SFE, ".025")

        PROPERTY_DOUBLE(cloudNumDensity, "Cloud number density")
        ATTRIBUTE_QUANTITY(cloudNumDensity, "numbervolumedensity")
        ATTRIBUTE_MIN_VALUE(cloudNumDensity, "[10")
        ATTRIBUTE_MAX_VALUE(cloudNumDensity, "2560]")
        ATTRIBUTE_DEFAULT_VALUE(cloudNumDensity, "320")

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

