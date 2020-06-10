/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef SINGLEWAVELENGTHSED_HPP
#define SINGLEWAVELENGTHSED_HPP

#include "SED.hpp"

////////////////////////////////////////////////////////////////////

/** The SingleWavelengthSED class implements a spectral energy distribution in the form of a
    Dirac-delta function. All photon packets are emitted at a single configurable wavelength,
    called the emission wavelength. The specific luminosity of the spectrum is undefined.
    Consequently,

    - This SED must be normalized through an integrated luminosity over a range that includes the
    emission wavelength rather than through a specific luminosity or a broadband.

    -The wavelength bias property of the corresponding source must be set to zero.

    The default value for the emission wavelength is the central wavelength of the hydrogen
    Lyman-alpha line. */
class SingleWavelengthSED : public SED
{
    ITEM_CONCRETE(SingleWavelengthSED, SED, "a single-wavelength SED in the form of a Dirac-delta function")

        PROPERTY_DOUBLE(wavelength, "the single emission wavelength")
        ATTRIBUTE_QUANTITY(wavelength, "wavelength")
        ATTRIBUTE_MIN_VALUE(wavelength, "1 Angstrom")
        ATTRIBUTE_MAX_VALUE(wavelength, "1 m")
        ATTRIBUTE_DEFAULT_VALUE(wavelength, "1215.67 Angstrom")

    ITEM_END()

    //======================== Other Functions =======================

public:
    /** This function returns the intrinsic wavelength range of the %SED. For this SED, the range
        just includes the emission wavelength. */
    Range intrinsicWavelengthRange() const override;

    /** This function returns the normalized specific luminosity \f$L_\lambda\f$ (i.e. radiative
        power per unit of wavelength) at the specified wavelength. For this SED, this value is
        undefined so the function throws a fatal error. */
    double specificLuminosity(double wavelength) const override;

    /** This function returns the normalized specific luminosity \f$L_\lambda\f$ (i.e. radiative
        power per unit of wavelength) at a number of wavelength points within the specified
        wavelength range. For this SED, these values are undefined so the function throws a fatal
        error. */
    void specificLuminosityArray(Array& lambdav, Array& pv, const Range& wavelengthRange) const override;

    /** This function returns the normalized integrated luminosity \f$L\f$ (i.e. radiative power)
        over the specified wavelength range. For this SED, the value is one if the specified range
        includes the emission wavelength and zero otherwise. */
    double integratedLuminosity(const Range& wavelengthRange) const override;

    /** This function draws a random wavelength from the normalized spectral energy distribution.
        For this SED, the emission wavelength is always returned. */
    double generateWavelength() const override;
};

////////////////////////////////////////////////////////////////////

#endif
