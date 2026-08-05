/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef MEANLISTDUSTMIX_HPP
#define MEANLISTDUSTMIX_HPP

#include "TabulatedDustMix.hpp"

////////////////////////////////////////////////////////////////////

/** A MeanListDustMix object represents a basic dust mix defined by optical properties fully
    specified inside the configuration file (i.e. without referring to an input file), describing a
    single representative grain and using the Henyey-Greenstein scattering mode. It is intended for
    use in cases where the properties are tabulated for just a few wavelengths (or even a single
    wavelength), but nothing keeps the user from specifying a long list.

    The numbers in the four user-configurable lists specify respectively the wavelength
    \f$\lambda\f$, the extinction mass coefficient \f$\kappa^\text{ext}_\lambda\f$, the scattering
    albedo \f$\varpi_\lambda\f$, and the scattering asymmetry parameter \f$g_\lambda\f$. The
    wavelengths must be listed in increasing or decreasing order. Property values outside of the
    tabulated wavelength range are clamped to the nearest border value. As a special-case
    consequence, if only a single wavelength is tabulated, the properties are considered to be
    constant for all wavelengths.

    Because a basic dust mix such as this one is usually used in isolation and the dust
    distribution is normalized to a given optical depth or total dust mass, the value of the
    extinction coefficient is essentially scale free. This class uses a fixed arbitrary value of
    the dust mass per hydrogen atom, \f$\mu=1.5\times 10^{-29} \text{kg}\,\text{H}^{-1}\f$ to set
    the absolute scale of the cross sections listed by some of the probes in relation to the mass
    coefficients. */
class MeanListDustMix : public TabulatedDustMix
{
    ITEM_CONCRETE(MeanListDustMix, TabulatedDustMix,
                  "a dust mix with mean properties specified inside the configuration file")

        PROPERTY_DOUBLE_LIST(wavelengths, "the wavelengths at which to specify the optical properties")
        ATTRIBUTE_QUANTITY(wavelengths, "wavelength")
        ATTRIBUTE_MIN_VALUE(wavelengths, "1 pm")
        ATTRIBUTE_MAX_VALUE(wavelengths, "1 m")

        PROPERTY_DOUBLE_LIST(extinctionCoefficients,
                             "the extinction mass coefficients at each of the given wavelengths")
        ATTRIBUTE_QUANTITY(extinctionCoefficients, "masscoefficient")
        ATTRIBUTE_MIN_VALUE(extinctionCoefficients, "[0")

        PROPERTY_DOUBLE_LIST(albedos, "the scattering albedos at each of the given wavelengths")
        ATTRIBUTE_MIN_VALUE(albedos, "[0")
        ATTRIBUTE_MAX_VALUE(albedos, "1]")

        PROPERTY_DOUBLE_LIST(asymmetryParameters,
                             "the scattering asymmetry parameters at each of the given wavelengths")
        ATTRIBUTE_MIN_VALUE(asymmetryParameters, "[-1")
        ATTRIBUTE_MAX_VALUE(asymmetryParameters, "1]")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function copies the user-configured properties into the specified arrays and returns a
        fixed, arbitrarily determined dust mass per hydrogen atom. */
    double getDustProperties(Array& lambdav, Array& kappaextv, Array& albedov, Array& asymmparv) const override;
};

////////////////////////////////////////////////////////////////////

#endif
