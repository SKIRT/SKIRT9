/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef STARBURST99SED_HPP
#define STARBURST99SED_HPP

#include "FamilySED.hpp"

////////////////////////////////////////////////////////////////////

/** Starburst99SED is a class that represents spectral energy distributions of simple stellar
    populations (SSPs), parameterized on metallicity and age according to the Starburst99 model
    using a Kroupa initial mass function. See the Starburst99SEDFamily class for more information.
    */
class Starburst99SED : public FamilySED
{
    ITEM_CONCRETE(Starburst99SED, FamilySED, "a Starburst99 simple stellar population SED")

        PROPERTY_DOUBLE(metallicity, "the metallicity of the SSP")
        ATTRIBUTE_MIN_VALUE(metallicity, "[0.0004")
        ATTRIBUTE_MAX_VALUE(metallicity, "0.05]")
        ATTRIBUTE_DEFAULT_VALUE(metallicity, "0.02")

        PROPERTY_DOUBLE(age, "the age of the SSP")
        ATTRIBUTE_QUANTITY(age, "time")
        ATTRIBUTE_MIN_VALUE(age, "[1e4 yr")
        ATTRIBUTE_MAX_VALUE(age, "1.53e10 yr]")
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
