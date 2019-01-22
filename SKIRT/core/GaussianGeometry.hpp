/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef GAUSSIANGEOMETRY_HPP
#define GAUSSIANGEOMETRY_HPP

#include "SpheGeometry.hpp"
#include "Array.hpp"

////////////////////////////////////////////////////////////////////

/** The GaussianGeometry class is a subclass of the SpheGeometry class, and describes spherical
    geometries characterized by a Gaussian density profile, \f[ \rho(r) = \rho_0\,\exp\left(
    -\frac{r^2}{2\sigma^2} \right) .\f] This geometry has one parameter, the radial dispersion
    \f$\sigma\f$. */
class GaussianGeometry : public SpheGeometry
{
    ITEM_CONCRETE(GaussianGeometry, SpheGeometry, "a Gaussian geometry")

    PROPERTY_DOUBLE(dispersion, "the scale length (dispersion) σ")
        ATTRIBUTE_QUANTITY(dispersion, "length")
        ATTRIBUTE_MIN_VALUE(dispersion, "]0")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function calculates some frequently used values. The central density \f$\rho_0\f$ is
        set by the normalization condition that the total mass equals one, which is straightforward
        for a Gaussian distribution, \f[ \rho_0 = \frac{1}{(2\pi)^{3/2}\,\sigma^3} .\f] This
        function also precalculates of a vector with the cumulative mass \f[ M(r) = 4\pi \int_0^r
        \rho(r')\, r'^2\, {\text{d}}r' \f] at a large number of radii. For the Gaussian
        distribution we find \f[ M(r) = \mathrm{erf}(t) - \frac{2}{\sqrt{\pi}} \,t \, \exp(-t^2)
        \quad\mathrm{with}\quad t = \frac{r}{\sqrt{2}\,\sigma} .\f] */
    void setupSelfBefore() override;

    //======================== Other Functions =======================

public:
    /** This function returns the density \f$\rho(r)\f$ at the radius \f$r\f$. It just implements
        the analytical formula. */
    double density(double r) const override;

    /** This function returns the radius \f$r\f$ of a random position drawn from a spherical
        Gaussian density distribution. Such a value can be generated by picking a uniform deviate
        \f${\cal{X}}\f$ and solving the equation \f[ {\cal{X}} = 4\pi \int_0^r \rho(r')\, r'^2\,
        {\text{d}}r' \f] for \f$r\f$, where \f$\rho(r)\f$ is the Gaussian radial density profile.
        This is done by interpolating from the precalculated table of the cumulative distribution.
        */
    double randomRadius() const override;

    /** This function returns the surface mass density along a radial line starting at the centre
        of the coordinate system, i.e. \f[ \Sigma_r = \int_0^\infty \rho(r)\,{\text{d}}r. \f] For a
        Gaussian geometry we easily find \f[ \Sigma_r = \frac{1}{4\pi\,\sigma^2}. \f] */
    double Sigmar() const override;

    //======================== Data Members ========================

private:
    // aliases to discoverable data members for ease of notation and backwards compatibility
    const double& _sigma{_dispersion};

    // data members initialized during setup
    double _rho0{0.};
    Array _rv;
    Array _Xv;
};

////////////////////////////////////////////////////////////////////

#endif
