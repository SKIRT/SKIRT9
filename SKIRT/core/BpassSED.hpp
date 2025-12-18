/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef BPASSSED_HPP
#define BPASSSED_HPP

#include "FamilySED.hpp"

////////////////////////////////////////////////////////////////////

/** BpassSED is a class that represents spectral energy distributions of simple stellar populations
    (SSPs), parameterized on metallicity and age according to the BPASS model that includes binary
    stellar systems and that assumes a Chabrier or Kroupa IMF with upper mass limits of 100 or 300 solar masses.
    See the BpassSEDFamily class for more information. */
class BpassSED : public FamilySED
{
    /** The enumeration type indicating the IMF to use */
    ENUM_DEF(IMF, Chabrier100, Chabrier300, Kroupa100, Kroupa300)
        ENUM_VAL(IMF, Chabrier100, "Chabrier IMF (0.1-100 Msun)")
        ENUM_VAL(IMF, Chabrier300, "Chabrier IMF (0.1-300 Msun)")
        ENUM_VAL(IMF, Kroupa100, "Kroupa IMF (0.1-100 Msun)")
        ENUM_VAL(IMF, Kroupa300, "Kroupa IMF (0.1-300 Msun)")
    ENUM_END()

    /** The enumeration type indicating the wavelength resolution and spectral range */
    ENUM_DEF(Resolution, Original, Downsampled)
        ENUM_VAL(Resolution, Original, "Original wavelength resolution and range")
        ENUM_VAL(Resolution, Downsampled, "Downsampled wavelength resolution and extended wavelength range")
    ENUM_END()

    ITEM_CONCRETE(BpassSED, FamilySED, "a BPASS single stellar population SED")

        PROPERTY_ENUM(imf, IMF, "the assumed initial mass function")
        ATTRIBUTE_DEFAULT_VALUE(imf, "Chabrier300")

        PROPERTY_ENUM(resolution, Resolution, "the wavelength resolution")
        ATTRIBUTE_DEFAULT_VALUE(resolution, "Original")

        PROPERTY_DOUBLE(metallicity, "the metallicity of the SSP")
        ATTRIBUTE_MIN_VALUE(metallicity, "[1e-5")
        ATTRIBUTE_MAX_VALUE(metallicity, "0.04]")
        ATTRIBUTE_DEFAULT_VALUE(metallicity, "0.02")

        PROPERTY_DOUBLE(age, "the age of the SSP")
        ATTRIBUTE_QUANTITY(age, "time")
        ATTRIBUTE_MIN_VALUE(age, "[1 Myr")
        ATTRIBUTE_MAX_VALUE(age, "100 Gyr]")
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
