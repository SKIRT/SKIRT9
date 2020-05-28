/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef FLATUNIVERSECOSMOLOGY_HPP
#define FLATUNIVERSECOSMOLOGY_HPP

#include "Cosmology.hpp"

////////////////////////////////////////////////////////////////////

/** An instance of the FlatUniverseCosmology class represents a cosmology where the simulated model
    resides at a given redshift \f$z\f$ in a standard spatially-flat \f$\Lambda\mathrm{CDM}\f$
    universe. The allowed redshift range is \f$0\le z\le 15\f$, i.e. including the epoch of
    reionization. If \f$z=0\f$, the behavior is the same as that of a LocalUniverseCosmology.

    In addition to the redshift, this cosmology has two parameters, namely the reduced Hubble
    constant \f$h=H_0/(100\,\mathrm{km}\,\mathrm{s}^{-1}\,\mathrm{Mpc}^{-1})\f$ and the matter
    density fraction \f$\Omega_\mathrm{m}\f$. The default values for these parameters are
    compatible with the latest observations and close to the values used by many recent
    cosmological simulations.

    For a spatially-flat cosmology, the line-of-sight and transverse comoving distances are
    identical. We further ignore the radiation density, which is justified for the allowed redshift
    range \f$z\le15\f$. With these assumptions, the comoving distance \f$d_\mathrm{M}(z)\f$
    corresponding to redshift $z$ can be obtained from the matter density fraction
    \f$\Omega_\mathrm{m}\f$ and the Hubble constant \f$H_0=h \times
    100\,\mathrm{km}\,\mathrm{s}^{-1}\,\mathrm{Mpc}^{-1}\f$ using \f[ d_\mathrm{M}(z) =
    \frac{c}{H_0} \int_0^z \frac{\mathrm{d}z'}{\sqrt{\Omega_\mathrm{m}(1+z')^3 +
    (1-\Omega_\mathrm{m})}} \f] where \f$c\f$ is the speed of light.

    The angular-diameter distance \f$d_\mathrm{A}(z)\f$ and the luminosity distance
    \f$d_\mathrm{L}(z)\f$ are then obtained from \f[ d_\mathrm{A}(z) = (1+z)^{-1} \,
    d_\mathrm{M}(z) \f] and \f[ d_\mathrm{L}(z) = (1+z) \, d_\mathrm{M}(z). \f]

    Finally, the relative expansion rate of the universe \f$H=\dot{a}/a\f$, usually called the
    Hubble parameter, is given by \f[ H(z) = H_0 \sqrt{\Omega_\mathrm{m}(1+z')^3 +
    (1-\Omega_\mathrm{m})}. \f]

    */
class FlatUniverseCosmology : public Cosmology
{
    ITEM_CONCRETE(FlatUniverseCosmology, Cosmology, "the model is at a given redshift in a flat universe")

        PROPERTY_DOUBLE(redshift, "the redshift z of the model coordinate frame")
        ATTRIBUTE_MIN_VALUE(redshift, "[0")
        ATTRIBUTE_MAX_VALUE(redshift, "15]")
        ATTRIBUTE_DEFAULT_VALUE(redshift, "1")
        ATTRIBUTE_INSERT(FlatUniverseCosmology, "redshift:NonZeroRedshift")

        PROPERTY_DOUBLE(reducedHubbleConstant, "the reduced Hubble constant h")
        ATTRIBUTE_MIN_VALUE(reducedHubbleConstant, "0.3")
        ATTRIBUTE_MAX_VALUE(reducedHubbleConstant, "1")
        ATTRIBUTE_DEFAULT_VALUE(reducedHubbleConstant, "0.675")
        ATTRIBUTE_RELEVANT_IF(reducedHubbleConstant, "redshift")

        PROPERTY_DOUBLE(matterDensityFraction, "the cosmological matter density fraction Ω_m")
        ATTRIBUTE_MIN_VALUE(matterDensityFraction, "0.1")
        ATTRIBUTE_MAX_VALUE(matterDensityFraction, "0.9")
        ATTRIBUTE_DEFAULT_VALUE(matterDensityFraction, "0.310")
        ATTRIBUTE_RELEVANT_IF(matterDensityFraction, "redshift")

    ITEM_END()

    //============== Functions overridden from base class ============

public:
    /** This function returns the redshift at which the model resides. */
    double modelRedshift() const override;

    /** This function returns the angular-diameter distance \f$d_\mathrm{A}(z)\f$, calculated from
        the user-configured cosmology parameters as described in the class header. The integral is
        evaluated numerically. We don't worry about performance because this function is called
        just once during simulation setup. */
    double angularDiameterDistance() const override;

    /** This function returns the luminosity distance. For this class, calculated from the
        user-configured cosmology parameters as described in the class header. The integral is
        evaluated numerically. We don't worry about performance because this function is called
        just once during simulation setup. */
    double luminosityDistance() const override;

    /** This function returns the relative expansion rate of the universe, \f$H=\dot{a}/a\f$, often
        called the Hubble parameter. The function returns the value in SI units, i.e.
        \f$\mathrm{s}^{-1}\f$. */
    double relativeExpansionRate() const override;
};

////////////////////////////////////////////////////////////////////

#endif
