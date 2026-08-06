/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef SERSICGEOMETRY_HPP
#define SERSICGEOMETRY_HPP

#include "SpheGeometry.hpp"
class SersicFunction;

////////////////////////////////////////////////////////////////////

/** The SersicGeometry class is a subclass of the SpheGeometry class, and describes spherically
    symmetric stellar geometries characterized by the density distribution \f[ \rho(r) = \rho_0\,
    {\cal{S}}_n \left(\frac{r}{r_{\text{eff}}}\right), \f] with \f${\cal{S}}_n(s)\f$ the Sersic
    function of order \f$n\f$ (see SersicFunction). It is defined in such a way that the projected
    surface brightness profile has the form \f[ I(r_p) = I_0 \exp \left[ -b_n\left(
    \frac{r_p}{r_{\text{eff}}} \right)^{1/n} \right]. \f] Two parameters characterize a
    SersicGeometry class object: the Sersic index \f$n\f$ and the effective radius
    \f$r_{\text{eff}}\f$. Internally, the SersicGeometry class has a SersicFunction object as one
    of its data members. See Sersic (1963) and Ciotti & Bertin (1999, A&A 352, 447–451). */
class SersicGeometry : public SpheGeometry
{
    ITEM_CONCRETE(SersicGeometry, SpheGeometry, "a Sérsic geometry")

        PROPERTY_DOUBLE(effectiveRadius, "the effective radius")
        ATTRIBUTE_QUANTITY(effectiveRadius, "length")
        ATTRIBUTE_MIN_VALUE(effectiveRadius, "]0")

        PROPERTY_DOUBLE(index, "the Sérsic index n")
        ATTRIBUTE_MIN_VALUE(index, "]0.5")
        ATTRIBUTE_MAX_VALUE(index, "10]")
        ATTRIBUTE_DEFAULT_VALUE(index, "1")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

public:
    /** The destructor; it deallocates the sersic function object created during setup. */
    ~SersicGeometry();

protected:
    /** This function creates a SersicFunction class object. The central density \f$\rho_0\f$ is
        set by the normalization condition that the total mass is equal to one. Since the Sersic
        function satisfies the normalization \f[ 4\pi \int_0^\infty {\cal{S}}_n(s)\, s^2\,
        {\text{d}}s = 1, \f] we easily find \f[ \rho_0 = \frac{1}{r_{\text{eff}}^3} \f] This
        function also caches the value of the dimensionless constant \f$b_n\f$ that appears in the
        definition of the Sersic profile. A suitable approximation is \f[ b_n = 2n -\frac{1}{3} +
        \frac{4}{405n} + \frac{46}{25515n^2} + \frac{131}{1148175n^3}. \f] */
    void setupSelfBefore() override;

    //======================== Other Functions =======================

public:
    /** This function returns the density \f$\rho(r)\f$ at a radius \f$r\f$. It just implements the
        analytical formula. */
    double density(double r) const override;

    /** This function returns the radius of a random position drawn from the Sersic density
        distribution. This is accomplished by generating a uniform deviate \f${\cal{X}}\f$, and
        solving the equation \f[ {\cal{X}} = M(r) = 4\pi \int_0^r \rho(r')\, r'{}^2\, {\text{d}}r'
        \f] for \f$r\f$. For the Sersic model, we use the SersicFunction::inversemass() function to
        solve this equation. */
    double randomRadius() const override;

    /** This function returns the radial surface density, i.e. the integration of
        the density along a line starting at the centre of the coordinate
        system, \f[ \Sigma_r = \int_0^\infty \rho(r)\,{\text{d}}r. \f] For the Sersic
        geometry, one finds \f[ \Sigma_r =
        \frac{1}{r_{\text{eff}}^2}\, \frac{b_n^{2n}}{ 2\pi\, \Gamma(2n+1)}. \f] */
    double Sigmar() const override;

    //======================== Data Members ========================

private:
    // aliases to discoverable data members for ease of notation and backwards compatibility
    const double& _reff{_effectiveRadius};
    const double& _n{_index};

    // data members initialized during setup
    double _rho0{0.};
    double _b{0.};
    SersicFunction* _sersicfunction{nullptr};
};

////////////////////////////////////////////////////////////////////

#endif
