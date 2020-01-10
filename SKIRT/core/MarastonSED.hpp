/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef MARASTONSED_HPP
#define MARASTONSED_HPP

#include "FamilySED.hpp"

////////////////////////////////////////////////////////////////////

/** MarastonSED is a class that represents spectral energy distributions of simple stellar
    populations (SSPs), parameterized on metallicity and age according to the Maraston 1998 model
    assuming a Kroupa or Salpeter initial mass function. See the MarastonSEDFamily class for more
    information. */
class MarastonSED : public FamilySED
{
    /** The enumeration type indicating the assumed initial mass function (IMF). */
    ENUM_DEF(IMF, Kroupa, Salpeter)
        ENUM_VAL(IMF, Kroupa, "Kroupa IMF")
        ENUM_VAL(IMF, Salpeter, "Salpeter IMF")
    ENUM_END()

    ITEM_CONCRETE(MarastonSED, FamilySED, "a Maraston simple stellar population SED")

        PROPERTY_ENUM(imf, IMF, "the assumed initial mass function")
        ATTRIBUTE_DEFAULT_VALUE(imf, "Kroupa")

        PROPERTY_DOUBLE(metallicity, "the metallicity of the SSP")
        ATTRIBUTE_MIN_VALUE(metallicity, "[0.0001")
        ATTRIBUTE_MAX_VALUE(metallicity, "0.09]")
        ATTRIBUTE_DEFAULT_VALUE(metallicity, "0.02")

        PROPERTY_DOUBLE(age, "the age of the SSP")
        ATTRIBUTE_QUANTITY(age, "time")
        ATTRIBUTE_MIN_VALUE(age, "[1000 yr")
        ATTRIBUTE_MAX_VALUE(age, "15 Gyr]")
        ATTRIBUTE_DEFAULT_VALUE(age, "5 Gyr")

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
