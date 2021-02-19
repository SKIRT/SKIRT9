/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef LYADOUBLEPEAKEDSED_HPP
#define LYADOUBLEPEAKEDSED_HPP

#include "ContSED.hpp"

////////////////////////////////////////////////////////////////////

/** The LyaDoublePeakedSED class implements a double-peaked %SED, centered on the Lyman-alpha line
    wavelength \f$\lambda_\alpha\f$, and corresponding to the spectrum that emerges from a static
    sphere of hydrogen gas surrounding a Lyman-alpha point source. This %SED can be used as a
    simple model to represent sub-grid gas embedded in the source.

    The spectrum can be expressed analytically (see Dijkstra et al. 2006a, ApJ, 649, 14-36). Using
    the photon velocity shift \f[ v = \frac{\lambda - \lambda_\alpha} {\lambda_\alpha} \,c \f] as
    the spectral variable and a velocity scale \f$s\f$ configured by the user, the normalized
    spectrum can be written as \f[ L_v(v) = \frac{3v^2} {2s^3\left[1+\cosh(v^3/s^3)\right]}. \f]
    The two peaks of this profile are situated at \f$v_* \approx \pm 1.06938 \,s\f$.

    The cumulative distribution is \f[ \int_{-\infty}^{v} L_v(v') \,\mathrm{d}v' = \frac{1}
    {1+\exp(-v^3/s^3)}, \f] which can be inverted so that we can obtain a sample from the
    distribution through \f[ v = s \left( \ln\frac{\mathcal{X}}{1-\mathcal{X}} \right)^{1/3}, \f]
    where \f$\mathcal{X}\f$ is a uniform deviate.

    To simplify some numerical operations, such as tabulating the %SED, we limit the intrinsic
    range of the %SED to \f$\pm 3.6s\f$ from the center. We have \f$L_v(3.6s)/L_v(v_*) <
    4\times10^{-19}\f$ so that values outside this range will never contribute to the results
    anyway.

    To simplify normalization of this %SED, the source wavelength range configured by the user must
    fully contain the intrinsic wavelength range defined by \f[ \lambda_\alpha \left(1 \pm
    \frac{3.6\,s}{c}\right) \f] */
class LyaDoublePeakedSED : public ContSED
{
    ITEM_CONCRETE(LyaDoublePeakedSED, ContSED, "a double-peaked spectrum around the central Lyman-alpha wavelength")
        ATTRIBUTE_TYPE_DISPLAYED_IF(LyaDoublePeakedSED, "Lya|Level3")

        PROPERTY_DOUBLE(scale, "the velocity scale")
        ATTRIBUTE_QUANTITY(scale, "velocity")
        ATTRIBUTE_MIN_VALUE(scale, "]0 km/s")
        ATTRIBUTE_MAX_VALUE(scale, "1000 km/s]")
        ATTRIBUTE_DEFAULT_VALUE(scale, "50 km/s")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function precalculates some values for later use. */
    void setupSelfBefore() override;

    //======================== Other Functions =======================

public:
    /** This function returns the intrinsic wavelength range of the %SED. For our double-peaked
        spectrum, we limit the range to \f$\pm 3.6s\f$ from the Lyman-alpha line center, as
        decribed in the class header. */
    Range intrinsicWavelengthRange() const override;

    /** This function returns the normalized specific luminosity \f$L_\lambda\f$ (i.e. radiative
        power per unit of wavelength) at the specified wavelength. */
    double specificLuminosity(double wavelength) const override;

    /** This function returns the normalized specific luminosity \f$L_\lambda\f$ (i.e. radiative
        power per unit of wavelength) at a number of wavelength points within the specified
        wavelength range. The intrinsic wavelength range \f$\pm 3.6s\f$ of the %SED is covered by
        720 wavelength points on a regular linear grid. */
    void specificLuminosityArray(Array& lambdav, Array& pv, const Range& wavelengthRange) const override;

    /** This function returns the normalized integrated luminosity \f$L\f$ (i.e. radiative power)
        over the specified wavelength range. */
    double integratedLuminosity(const Range& wavelengthRange) const override;

    /** This function draws a random wavelength from the normalized spectral energy distribution.
        */
    double generateWavelength() const override;

    //======================== Data Members ========================

private:
    double _wavelengthCenter{0};  // the center of the double-peaked spectrum in wavelength units
    double _wavelengthScale{0};   // the scale of the double-peaked spectrum in wavelength units
    Range _wavelengthRange;       // the assumed intrinsic range in wavelength units
};

////////////////////////////////////////////////////////////////////

#endif
