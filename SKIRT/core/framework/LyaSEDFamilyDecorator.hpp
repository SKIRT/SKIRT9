/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef LYASEDFAMILYDECORATOR_HPP
#define LYASEDFAMILYDECORATOR_HPP

#include "ContSED.hpp"
#include "SEDFamily.hpp"

////////////////////////////////////////////////////////////////////

/** The LyaSEDFamilyDecorator class implements a \em decorator that adjusts another, arbitrary %SED
    family by converting a fraction of the ionizing part of the %SEDs to Lyman-alpha emission.

    There are three configuration options: the %SED family to be adjusted, called the 'original
    %SED family'; the %SED representing the Lyman-alpha emission, called the 'Lyman-alpha %SED';
    and the fraction of the ionizing radiation to be converted. The template parameters passed to
    the decorator are passed on to the original %SED family without change. In the 'decorated %SED'
    templates produced by the decorator, the specified fraction of the luminosity in any original
    %SED short of \f$\lambda_\mathrm{ion} =911.75\,\text{\AA}\f$ is replaced by emission following
    the specified Lyman-alpha %SED. The remaining fraction of ionizing radiation and all
    non-ionizing radiation in the original %SED is emitted as usual.

    The source wavelength range must include the ionizing portion of the spectrum to be considered,
    even if all ionizing radiation is being converted to Lyman-alpha emission. */
class LyaSEDFamilyDecorator : public SEDFamily
{
    ITEM_CONCRETE(LyaSEDFamilyDecorator, SEDFamily,
                  "an SED family decorator replacing ionizing radiation with Lyman-alpha line emission")
        ATTRIBUTE_TYPE_DISPLAYED_IF(LyaSEDFamilyDecorator, "Lya|Level3")

        PROPERTY_ITEM(sedFamilyOriginal, SEDFamily, "the original SED family being adjusted")
        ATTRIBUTE_DEFAULT_VALUE(sedFamilyOriginal, "BlackBodySEDFamily")

        PROPERTY_ITEM(sedLymanAlpha, ContSED, "the SED representing the Lyman-alpha line emission")
        ATTRIBUTE_DEFAULT_VALUE(sedLymanAlpha, "LyaGaussianSED")

        PROPERTY_DOUBLE(conversionFraction, "the fraction of ionizing radiation replaced by Lyman-alpha emission")
        ATTRIBUTE_MIN_VALUE(conversionFraction, "[0")
        ATTRIBUTE_MAX_VALUE(conversionFraction, "1]")
        ATTRIBUTE_DEFAULT_VALUE(conversionFraction, "0.5")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function verifies that the source wavelength range includes ionizing radiation and
        caches the ionizing wavelength range. */
    void setupSelfAfter() override;

    //======================== Other Functions =======================

public:
    /** This function passes on the number and type of parameters used by the original %SED family
        as a list of SnapshotParameter objects. Each of these objects specifies unit information
        and a human-readable descripton for the parameter. */
    vector<SnapshotParameter> parameterInfo() const override;

    /** This function returns the intrinsic wavelength range of the %SED family. For this
        decorator, it returns the union of the intrinsic ranges of the original %SED family and the
        Lyman-alpha %SED configured by the user, even if the configured conversion fraction is
        equal to one. */
    Range intrinsicWavelengthRange() const override;

    /** This function returns the specific luminosity \f$L_\lambda\f$ (i.e. radiative power per
        unit of wavelength) for the decorated %SED with the specified parameters at the specified
        wavelength, or zero if the wavelength is outside of the %SED's intrinsic wavelength range.
        The number and type of parameters must match the information returned by the
        parameterInfo() function (which match those of the original %SED family); if not the
        behavior is undefined. */
    double specificLuminosity(double wavelength, const Array& parameters) const override;

    /** This function constructs both the normalized probability density function (pdf) and the
        corresponding normalized cumulative distribution function (cdf) for the decorated %SED with
        the specified parameters over the specified wavelength range. The function returns the
        normalization factor. The number and type of parameters must match the information returned
        by the parameterInfo() function (which match those of the original %SED family); if not the
        behavior is undefined. */
    double cdf(Array& lambdav, Array& pv, Array& Pv, const Range& wavelengthRange,
               const Array& parameters) const override;

    //======================== Data Members ========================

private:
    Range _ionizingRange;  // the wavelength range representing ionizing radiation
};

////////////////////////////////////////////////////////////////////

#endif
