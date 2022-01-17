/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef LINELUMINOSITYNORMALIZATION_HPP
#define LINELUMINOSITYNORMALIZATION_HPP

#include "LuminosityNormalization.hpp"

////////////////////////////////////////////////////////////////////

/** A LineLuminosityNormalization instance sets the normalization of a primary source by specifying
    the luminosity (radiative power) for a given emission line. In practice, the normalization
    occurs by integrating over a very narrow wavelength range centered on a given wavelength value.
    This class is intended for use with sources of line emission, i.e. sources equipped with one of
    the LineSED subclasses. */
class LineLuminosityNormalization : public LuminosityNormalization
{
    ITEM_CONCRETE(LineLuminosityNormalization, LuminosityNormalization,
                  "source normalization through the luminosity for a given emission line")
        ATTRIBUTE_TYPE_ALLOWED_IF(LineLuminosityNormalization, "LineSED")
        ATTRIBUTE_TYPE_DISPLAYED_IF(LineLuminosityNormalization, "Level2")

        PROPERTY_DOUBLE(wavelength, "the central wavelength of the emission line")
        ATTRIBUTE_QUANTITY(wavelength, "wavelength")
        ATTRIBUTE_MIN_VALUE(wavelength, "1 pm")
        ATTRIBUTE_MAX_VALUE(wavelength, "1 m")
        ATTRIBUTE_DEFAULT_VALUE(wavelength, "1215.67 Angstrom")

        PROPERTY_DOUBLE(luminosity, "the luminosity for the emission line")
        ATTRIBUTE_QUANTITY(luminosity, "bolluminosity")
        ATTRIBUTE_MIN_VALUE(luminosity, "]0")

    ITEM_END()

    //======================== Other Functions =======================

public:
    /** This function returns the luminosity of a source with the specified %SED, limited to the
        source's wavelength range. Given that the luminosity of the spectrum described by the SED
        object over the source wavelength range is normalized to unity, the requested luminosity is
        obtained by dividing the user-configured luminosity by the normalized luminosity for the
        line in the specified %SED. */
    double luminosityForSED(SED* sed) const override;
};

////////////////////////////////////////////////////////////////////

#endif
