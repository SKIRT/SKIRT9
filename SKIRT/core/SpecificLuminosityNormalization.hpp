/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef SPECTRALLUMINOSITYSTELLARCOMPNORMALIZATION_HPP
#define SPECTRALLUMINOSITYSTELLARCOMPNORMALIZATION_HPP

#include "LuminosityNormalization.hpp"

////////////////////////////////////////////////////////////////////

/** SpectralLuminosityStellarCompNormalization is a class that sets the normalization of a stellar
    component by defining the spectral luminosity (radiative power per wavelength) at a certain
    wavelength. */
class SpecificLuminosityNormalization : public LuminosityNormalization
{
    ITEM_CONCRETE(SpecificLuminosityNormalization, LuminosityNormalization,
                  "source normalization through the specific luminosity at a given wavelength")

    PROPERTY_DOUBLE(wavelength, "the wavelength at which to specify the specific luminosity")
        ATTRIBUTE_QUANTITY(wavelength, "wavelength")
        ATTRIBUTE_MIN_VALUE(wavelength, "1 A")
        ATTRIBUTE_MAX_VALUE(wavelength, "1 m")

    PROPERTY_DOUBLE(specificLuminosity, "the specific luminosity for the source at the specified wavelength")
        ATTRIBUTE_QUANTITY(specificLuminosity, "wavelengthmonluminosity")
        ATTRIBUTE_MIN_VALUE(specificLuminosity, "]0")

    ITEM_END()

    //======================== Other Functions =======================

public:
    /** This function returns the bolometric luminosity of a source with the specified %SED. For
        the present type of normalization, the bolometric luminosity is obtained by dividing the
        user-configured specific luminosity by the normalized specific luminosity in the %SED at
        the same wavelength. This follows from the fact that the bolometric luminosity of the
        spectrum described by the SED object is one by definition. */
    double luminosity(SED* sed) const override;
};

////////////////////////////////////////////////////////////////////

#endif
