/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef LYASEDDECORATOR_HPP
#define LYASEDDECORATOR_HPP

#include "ContSED.hpp"

////////////////////////////////////////////////////////////////////

/** The LyaSEDDecorator class implements a \em decorator that adjusts another, arbitrary %SED by
    converting a fraction of the ionizing part of that %SED to Lyman-alpha emission.

    There are three configuration options: the %SED to be adjusted, called the 'original %SED'; the
    %SED representing the Lyman-alpha emission, called the 'Lyman-alpha %SED'; and the fraction of
    the ionizing radiation to be converted. In the 'decorated %SED' produced by this class, the
    specified fraction of the luminosity in the original %SED short of \f$\lambda_\mathrm{ion}
    =911.75\,\text{\AA}\f$ is replaced by emission following the specified Lyman-alpha %SED. The
    remaining fraction of ionizing radiation and all non-ionizing radiation in the original %SED is
    emitted as usual.

    The luminosity normalization configured by the user applies to the decorated %SED. The source
    wavelength range must include the ionizing portion of the spectrum to be considered, even if
    all ionizing radiation is being converted to Lyman-alpha emission. */
class LyaSEDDecorator : public ContSED
{
    ITEM_CONCRETE(LyaSEDDecorator, ContSED,
                  "an SED decorator replacing ionizing radiation with Lyman-alpha line emission")
        ATTRIBUTE_TYPE_DISPLAYED_IF(LyaSEDDecorator, "Lya|Level3")

        PROPERTY_ITEM(sedOriginal, ContSED, "the original SED being adjusted")
        ATTRIBUTE_DEFAULT_VALUE(sedOriginal, "BlackBodySED")

        PROPERTY_ITEM(sedLymanAlpha, ContSED, "the SED representing the Lyman-alpha line emission")
        ATTRIBUTE_DEFAULT_VALUE(sedLymanAlpha, "LyaGaussianSED")

        PROPERTY_DOUBLE(conversionFraction, "the fraction of ionizing radiation replaced by Lyman-alpha emission")
        ATTRIBUTE_MIN_VALUE(conversionFraction, "[0")
        ATTRIBUTE_MAX_VALUE(conversionFraction, "1]")
        ATTRIBUTE_DEFAULT_VALUE(conversionFraction, "0.5")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function precalculates some values for later use. */
    void setupSelfAfter() override;

    //======================== Other Functions =======================

public:
    /** This function returns the intrinsic wavelength range of the %SED. For this decorator, it
        returns the union of the intrinsic ranges of the original and Lyman-alpha SEDs configured by
        the user, even if the configured conversion fraction is equal to one. */
    Range intrinsicWavelengthRange() const override;

    /** This function returns the normalized specific luminosity \f$L_\lambda\f$ (i.e. radiative
        power per unit of wavelength) at the specified wavelength. It determines the appropriate
        value from the user-configured SEDs as described in the class header. */
    double specificLuminosity(double wavelength) const override;

    /** This function returns the normalized specific luminosity \f$L_\lambda\f$ (i.e. radiative
        power per unit of wavelength) at a number of wavelength points within the specified
        wavelength range. It determines the appropriate values from the user-configured SEDs as
        described in the class header. */
    void specificLuminosityArray(Array& lambdav, Array& pv, const Range& wavelengthRange) const override;

    /** This function returns the normalized integrated luminosity \f$L\f$ (i.e. radiative power)
        over the specified wavelength range. It determines the appropriate value from the
        user-configured SEDs as described in the class header. */
    double integratedLuminosity(const Range& wavelengthRange) const override;

    /** This function draws a random wavelength from the normalized spectral energy distribution.
        It generates a wavelength from one of the user-configured SEDs as described in the class
        header. */
    double generateWavelength() const override;

    //======================== Data Members ========================

private:
    double _ionizingFraction{0};  // the fraction of the original luminosity representing ionizing radiation
};

////////////////////////////////////////////////////////////////////

#endif
