/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef MEANFILEDUSTMIX_HPP
#define MEANFILEDUSTMIX_HPP

#include "TabulatedDustMix.hpp"

////////////////////////////////////////////////////////////////////

/** A MeanFileDustMix object represents a basic dust mix defined by optical properties loaded from
    a file, describing a single representative grain and using the Henyey-Greenstein scattering
    mode.

    The floating point numbers in the first four columns of the text file specify respectively the
    wavelength \f$\lambda\f$, the extinction mass coefficient \f$\kappa^\text{ext}_\lambda\f$, the
    scattering albedo \f$\varpi_\lambda\f$, and the scattering asymmetry parameter \f$g_\lambda\f$.
    Any additional columns in the file are ignored.

    The wavelengths are by default given in micron (the units can be overridden by column header
    info in the file). The extinction mass coefficients are by default given in
    \f$\text{m}^2\,\text{kg}^{-2}\f$ (again, this can be overridden by column header info in the
    file). The wavelengths must be listed in increasing or decreasing order. Property values
    outside of the tabulated wavelength range are clamped to the nearest border value. As a
    special-case consequence, if only a single wavelength is tabulated, the properties are
    considered to be constant for all wavelengths.

    Because a basic dust mix such as this one is usually used in isolation and the dust
    distribution is normalized to a given optical depth or total dust mass, the value of the
    extinction coefficient is essentially scale free. This class uses a fixed arbitrary value of
    the dust mass per hydrogen atom, \f$\mu=1.5\times 10^{-29} \text{kg}\,\text{H}^{-1}\f$ to set
    the absolute scale of the cross sections listed by some of the probes in relation to the mass
    coefficients. */
class MeanFileDustMix : public TabulatedDustMix
{
    ITEM_CONCRETE(MeanFileDustMix, TabulatedDustMix, "a dust mix with mean properties loaded from a text file")

        PROPERTY_STRING(filename, "the name of the file with the optical properties for the dust mix")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function loads the input file into the specified arrays and returns a fixed,
        arbitrarily determined dust mass per hydrogen atom. */
    double getDustProperties(Array& lambdav, Array& kappaextv, Array& albedov, Array& asymmparv) const override;
};

////////////////////////////////////////////////////////////////////

#endif
