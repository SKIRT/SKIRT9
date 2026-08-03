/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef LYAGAUSSIANSED_HPP
#define LYAGAUSSIANSED_HPP

#include "ContSED.hpp"

////////////////////////////////////////////////////////////////////

/** The LyaGaussianSED class implements a Gaussian spectrum around the central Lyman-alpha
    wavelength \f$\lambda_\alpha\f$, reflecting the thermal sub-grid motion in the source. Using
    the photon velocity shift \f[ v = \frac{\lambda - \lambda_\alpha} {\lambda_\alpha} \,c \f] as
    the spectral variable, and a dispersion \f$s\f$ configured by the user in velocity units, the
    normalized Gaussian spectrum can be written as \f[ L_v(v) = \frac{1}{s\,\sqrt{2\pi}}
    \,\exp\left( -\frac{v^2}{2s^2} \right). \f]

    The thermal velocity corresponding to the dispersion \f$s\f$ is given by \f[
    v_\mathrm{th}=\sqrt{2}\,s \f] and the full width at half maximum of the Gaussian line profile
    is given by \f[ \text{FWHM} = 2\sqrt{2\ln2}\,s \approx 2.35482\,s. \f]

    This %SED can be implemented analytically and it can sampled using the standard procedure for a
    Gaussian deviate. To simplify some numerical operations, such as tabulating the %SED, we limit
    the intrinsic range of the %SED to \f$\pm 9s\f$ from the center. We have \f$L_v(9s)/L_v(0) <
    3\times10^{-18}\f$ so that values outside this range will never contribute to the results
    anyway.

    To simplify normalization of this %SED, the source wavelength range configured by the user must
    fully contain the intrinsic wavelength range defined by \f[ \lambda_\alpha \left(1 \pm
    \frac{9\,s}{c}\right) \f] */
class LyaGaussianSED : public ContSED
{
    ITEM_CONCRETE(LyaGaussianSED, ContSED, "a Gaussian spectrum around the central Lyman-alpha wavelength")
        ATTRIBUTE_TYPE_DISPLAYED_IF(LyaGaussianSED, "Lya|Level3")

        PROPERTY_DOUBLE(dispersion, "the Gaussian velocity dispersion")
        ATTRIBUTE_QUANTITY(dispersion, "velocity")
        ATTRIBUTE_MIN_VALUE(dispersion, "]0 km/s")
        ATTRIBUTE_MAX_VALUE(dispersion, "1000 km/s]")
        ATTRIBUTE_DEFAULT_VALUE(dispersion, "50 km/s")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function precalculates some values for later use. */
    void setupSelfBefore() override;

    //======================== Other Functions =======================

public:
    /** This function returns the intrinsic wavelength range of the %SED. For our Gaussian profile,
        we limit the range to \f$\pm 9s\f$ from the Lyman-alpha line center, as decribed in the
        class header. */
    Range intrinsicWavelengthRange() const override;

    /** This function returns the normalized specific luminosity \f$L_\lambda\f$ (i.e. radiative
        power per unit of wavelength) at the specified wavelength. */
    double specificLuminosity(double wavelength) const override;

    /** This function returns the normalized specific luminosity \f$L_\lambda\f$ (i.e. radiative
        power per unit of wavelength) at a number of wavelength points within the specified
        wavelength range. The intrinsic wavelength range \f$\pm 9s\f$ of the %SED is covered by
        1800 wavelength points on a regular linear grid. */
    void specificLuminosityArray(Array& lambdav, Array& pv, const Range& wavelengthRange) const override;

    /** This function returns the normalized integrated luminosity \f$L\f$ (i.e. radiative power)
        over the specified wavelength range. */
    double integratedLuminosity(const Range& wavelengthRange) const override;

    /** This function draws a random wavelength from the normalized spectral energy distribution.
        */
    double generateWavelength() const override;

    //======================== Data Members ========================

private:
    double _wavelengthCenter{0};      // the center of the Gaussian distribution in wavelength units
    double _wavelengthDispersion{0};  // the dispersion of the Gaussian distribution in wavelength units
    Range _wavelengthRange;           // the assumed intrinsic range in wavelength units
};

////////////////////////////////////////////////////////////////////

#endif
