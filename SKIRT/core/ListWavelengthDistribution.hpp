/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef LISTWAVELENGTHDISTRIBUTION_HPP
#define LISTWAVELENGTHDISTRIBUTION_HPP

#include "TabulatedWavelengthDistribution.hpp"

////////////////////////////////////////////////////////////////////

/** A ListWavelengthDistribution object represents a wavelength probability distribution that is
    fully specified inside the configuration file (i.e. without referring to an input file) by
    tabulated wavelength/probability pairs. The probability distribution function is defined
    segment-wise by the tabulated values, using logarithmic interpolation. This class is intended
    for use in cases where there are just a few wavelength/probability pairs, but nothing keeps the
    user from specifying a long list. The wavelengths must be listed in increasing or decreasing
    order. The probability outside the range indicated by the first and the last wavelength in the
    list is considered to be zero.

    The probability values are in fact given in luminosity units. The default is to use
    per-wavelength units, but the user can opt to use per-frequency or neutral units. Other than
    this, the scaling of the values is arbitrary because the distribution will be normalized after
    being loaded. However, the input procedure still insists on knowing the precise units. */
class ListWavelengthDistribution : public TabulatedWavelengthDistribution
{
    /** The enumeration type indicating the specific luminosity unit style, e.g. whether to use
        specific luminosity per unit of wavelength, frequency or energy. */
    ENUM_DEF(UnitStyle, neutralmonluminosity, wavelengthmonluminosity, frequencymonluminosity, energymonluminosity)
        ENUM_VAL(UnitStyle, neutralmonluminosity, "neutral: λ L_λ = ν L_ν")
        ENUM_VAL(UnitStyle, wavelengthmonluminosity, "per unit of wavelength: L_λ")
        ENUM_VAL(UnitStyle, frequencymonluminosity, "per unit of frequency: L_ν")
        ENUM_VAL(UnitStyle, energymonluminosity, "counts per unit of energy: L_E")
    ENUM_END()

    ITEM_CONCRETE(ListWavelengthDistribution, TabulatedWavelengthDistribution,
                  "a wavelength probability distribution specified inside the configuration file")

        PROPERTY_DOUBLE_LIST(wavelengths, "the wavelengths at which to specify the probability")
        ATTRIBUTE_QUANTITY(wavelengths, "wavelength")
        ATTRIBUTE_MIN_VALUE(wavelengths, "1 pm")
        ATTRIBUTE_MAX_VALUE(wavelengths, "1 m")

        PROPERTY_ENUM(unitStyle, UnitStyle, "the probability unit style")
        ATTRIBUTE_DEFAULT_VALUE(unitStyle, "wavelengthmonluminosity")

        PROPERTY_DOUBLE_LIST(probabilities, "the probabilities at each of the given wavelengths")
        ATTRIBUTE_QUANTITY(probabilities, "@unitStyle")
        ATTRIBUTE_MIN_VALUE(probabilities, "]0")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function loads the configured data into the specified arrays. */
    void getWavelengthsAndProbabilities(Array& lambdav, Array& pv) const override;
};

////////////////////////////////////////////////////////////////////

#endif
