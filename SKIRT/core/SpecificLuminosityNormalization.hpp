/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef SPECTRALLUMINOSITYSTELLARCOMPNORMALIZATION_HPP
#define SPECTRALLUMINOSITYSTELLARCOMPNORMALIZATION_HPP

#include "LuminosityNormalization.hpp"

////////////////////////////////////////////////////////////////////

/** A SpecificLuminosityNormalization instance sets the normalization of a primary source by
    specifying the specific luminosity (radiative power per units of wavelength or frequency) at a
    certain wavelength. */
class SpecificLuminosityNormalization : public LuminosityNormalization
{
    ITEM_CONCRETE(SpecificLuminosityNormalization, LuminosityNormalization,
                  "source normalization through the specific luminosity at a given wavelength")

    PROPERTY_DOUBLE(wavelength, "the wavelength at which to provide the specific luminosity")
        ATTRIBUTE_QUANTITY(wavelength, "wavelength")
        ATTRIBUTE_MIN_VALUE(wavelength, "1 A")
        ATTRIBUTE_MAX_VALUE(wavelength, "1 m")

    PROPERTY_DOUBLE(specificLuminosity, "the specific luminosity at the given wavelength")
        ATTRIBUTE_QUANTITY(specificLuminosity, "wavelengthmonluminosity")
        ATTRIBUTE_MIN_VALUE(specificLuminosity, "]0")

    ITEM_END()

    //======================== Other Functions =======================

public:
    /** This function returns the luminosity of a source with the specified %SED, limited to the
        source's wavelength range. Given that the luminosity of the spectrum described by the SED
        object over the source wavelength range is normalized to unity, the requested luminosity is
        obtained by dividing the user-configured specific luminosity by the normalized specific
        luminosity in the %SED at the same wavelength. */
    double luminosity(SED* sed) const override;
};

////////////////////////////////////////////////////////////////////

#endif
