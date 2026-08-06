/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef CastelliKuruczSED_HPP
#define CastelliKuruczSED_HPP

#include "FamilySED.hpp"

////////////////////////////////////////////////////////////////////

/** CastelliKuruczSED is a class that represents spectral energy distributions of stellar
    atmospheres, parameterized on metallicity, effective temperature, and gravity, according to the
    Castelli-Kurucz 2003 models. See the CastelliKuruczSEDFamily class for more information. */
class CastelliKuruczSED : public FamilySED
{
    ITEM_CONCRETE(CastelliKuruczSED, FamilySED, "a Castelli-Kurucz stellar atmosphere SED")

        PROPERTY_DOUBLE(metallicity, "the stellar metallicity")
        ATTRIBUTE_MIN_VALUE(metallicity, "[0.00006")
        ATTRIBUTE_MAX_VALUE(metallicity, "0.06]")
        ATTRIBUTE_DEFAULT_VALUE(metallicity, "0.02")

        PROPERTY_DOUBLE(temperature, "the effective surface temperature")
        ATTRIBUTE_QUANTITY(temperature, "temperature")
        ATTRIBUTE_MIN_VALUE(temperature, "]3500 K")
        ATTRIBUTE_MAX_VALUE(temperature, "50000 K]")
        ATTRIBUTE_DEFAULT_VALUE(temperature, "5000 K")

        PROPERTY_DOUBLE(gravity, "the surface gravity")
        ATTRIBUTE_QUANTITY(gravity, "acceleration")
        ATTRIBUTE_MIN_VALUE(gravity, "[0.01 m/s2")
        ATTRIBUTE_MAX_VALUE(gravity, "1000 m/s2]")
        ATTRIBUTE_DEFAULT_VALUE(gravity, "275 m/s2")

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
