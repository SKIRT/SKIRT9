/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef INTEGRATEDLUMINOSITYNORMALIZATION_HPP
#define INTEGRATEDLUMINOSITYNORMALIZATION_HPP

#include "LuminosityNormalization.hpp"

////////////////////////////////////////////////////////////////////

/** An IntegratedLuminosityNormalization instance sets the normalization of a primary source by
    specifying the luminosity (radiative power) integrated over a certain wavelength range. */
class IntegratedLuminosityNormalization : public LuminosityNormalization
{
    /** The enumeration type selecting a wavelength range: the range given for the source system,
        the bolometric range of the SED (i.e. including all wavelengths), or a custom range
        specified as properties of this object. */
    ENUM_DEF(WavelengthRange, Source, All, Custom)
        ENUM_VAL(WavelengthRange, Source, "the wavelength range of the primary sources")
        ENUM_VAL(WavelengthRange, All, "all wavelengths (i.e. over the full SED)")
        ENUM_VAL(WavelengthRange, Custom, "a custom wavelength range specified here")
    ENUM_END()

    ITEM_CONCRETE(IntegratedLuminosityNormalization, LuminosityNormalization,
                  "source normalization through the integrated luminosity for a given wavelength range")

        PROPERTY_ENUM(wavelengthRange, WavelengthRange,
                      "the wavelength range for which to provide the integrated luminosity")
        ATTRIBUTE_DEFAULT_VALUE(wavelengthRange, "Source")

        PROPERTY_DOUBLE(minWavelength, "the shortest wavelength of the integration range")
        ATTRIBUTE_QUANTITY(minWavelength, "wavelength")
        ATTRIBUTE_MIN_VALUE(minWavelength, "1 pm")
        ATTRIBUTE_MAX_VALUE(minWavelength, "1 m")
        ATTRIBUTE_DEFAULT_VALUE(minWavelength, "0.09 micron")
        ATTRIBUTE_RELEVANT_IF(minWavelength, "wavelengthRangeCustom")

        PROPERTY_DOUBLE(maxWavelength, "the longest wavelength of the integration range")
        ATTRIBUTE_QUANTITY(maxWavelength, "wavelength")
        ATTRIBUTE_MIN_VALUE(maxWavelength, "1 pm")
        ATTRIBUTE_MAX_VALUE(maxWavelength, "1 m")
        ATTRIBUTE_DEFAULT_VALUE(maxWavelength, "100 micron")
        ATTRIBUTE_RELEVANT_IF(maxWavelength, "wavelengthRangeCustom")

        PROPERTY_DOUBLE(integratedLuminosity, "the integrated luminosity for the given wavelength range")
        ATTRIBUTE_QUANTITY(integratedLuminosity, "bolluminosity")
        ATTRIBUTE_MIN_VALUE(integratedLuminosity, "]0")

    ITEM_END()

    //======================== Other Functions =======================

public:
    /** This function returns the luminosity of a source with the specified %SED, limited to the
        source's wavelength range. Given that the luminosity of the spectrum described by the SED
        object over the source wavelength range is normalized to unity, the requested luminosity is
        obtained by dividing the user-configured integrated luminosity by the normalized integrated
        luminosity for the %SED over the same wavelength range. */
    double luminosityForSED(SED* sed) const override;
};

////////////////////////////////////////////////////////////////////

#endif
