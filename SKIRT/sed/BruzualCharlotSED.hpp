/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef BRUZUALCHARLOTSED_HPP
#define BRUZUALCHARLOTSED_HPP

#include "FamilySED.hpp"

////////////////////////////////////////////////////////////////////

/** BruzualCharlotSED is a class that represents spectral energy distributions of simple stellar
    populations (SSPs), parameterized on metallicity and age according to the Bruzual & Charlot
    2003 model assuming a Chabrier or Salpeter initial mass function. It always uses the high
    wavelength resolution version. See the BruzualCharlotSEDFamily class for more information. */
class BruzualCharlotSED : public FamilySED
{
    /** The enumeration type indicating the assumed initial mass function (IMF). */
    ENUM_DEF(IMF, Chabrier, Salpeter)
        ENUM_VAL(IMF, Chabrier, "Chabrier IMF")
        ENUM_VAL(IMF, Salpeter, "Salpeter IMF")
    ENUM_END()

    ITEM_CONCRETE(BruzualCharlotSED, FamilySED, "a Bruzual-Charlot simple stellar population SED")

        PROPERTY_ENUM(imf, IMF, "the assumed initial mass function")
        ATTRIBUTE_DEFAULT_VALUE(imf, "Chabrier")

        PROPERTY_DOUBLE(metallicity, "the metallicity of the SSP")
        ATTRIBUTE_MIN_VALUE(metallicity, "[0.0001")
        ATTRIBUTE_MAX_VALUE(metallicity, "0.05]")
        ATTRIBUTE_DEFAULT_VALUE(metallicity, "0.02")

        PROPERTY_DOUBLE(age, "the age of the SSP")
        ATTRIBUTE_QUANTITY(age, "time")
        ATTRIBUTE_MIN_VALUE(age, "[0 Gyr")
        ATTRIBUTE_MAX_VALUE(age, "20 Gyr]")
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
