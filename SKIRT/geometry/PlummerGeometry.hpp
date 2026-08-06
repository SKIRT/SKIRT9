/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef PLUMMERGEOMETRY_HPP
#define PLUMMERGEOMETRY_HPP

#include "SpheGeometry.hpp"

////////////////////////////////////////////////////////////////////

/** The PlummerGeometry class is a subclass of the SpheGeometry class, and describes spherically
    symmetric geometries characterized by a Plummer density profile, \f[ \rho(r) = \rho_0\,
    \left(1+\frac{r^2}{c^2}\right)^{-5/2}. \f] The only free parameter is the scale length \f$c\f$.
    See Plummer (1911, MNRAS, 71, 460–470) and Dejonghe (1987, MNRAS, 224, 13–39). */
class PlummerGeometry : public SpheGeometry
{
    ITEM_CONCRETE(PlummerGeometry, SpheGeometry, "a Plummer geometry")

        PROPERTY_DOUBLE(scaleLength, "the scale length")
        ATTRIBUTE_QUANTITY(scaleLength, "length")
        ATTRIBUTE_MIN_VALUE(scaleLength, "]0")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function calculates some frequently used values. The central density
        \f$\rho_0\f$ is set by the normalization condition that the total mass equals one. For the
        Plummer model we find \f[ \rho_0 = \frac{3}{4\pi\,c^3}.\f] */
    void setupSelfBefore() override;

    //======================== Other Functions =======================

public:
    /** This function returns the density \f$\rho(r)\f$ at a radius \f$r\f$. It just implements the
        analytical formula. */
    double density(double r) const override;

    /** This function returns the radius of a random position drawn from the Plummer density
        distribution. This is accomplished by generating a uniform deviate \f${\cal{X}}\f$, and
        solving the equation \f[ {\cal{X}} = M(r) = 4\pi \int_0^r \rho(r')\, r'{}^2\, {\text{d}}r'
        \f] for \f$r\f$. For the Plummer model, we obtain the simple expression \f[ r = c\,
        \frac{{\cal{X}}^{1/3}}{\sqrt{1-{\cal{X}}^{2/3}}}. \f] */
    double randomRadius() const override;

    /** This function returns the radial surface density, i.e. the integration of
        the density along a line starting at the centre of the coordinate
        system, \f[ \Sigma_r = \int_0^\infty \rho(r)\,{\text{d}}r. \f] For the Plummer geometry,
        one obtains \f[ \Sigma_r = \frac{1}{2\pi\,c^2}. \f] */
    double Sigmar() const override;

    //======================== Data Members ========================

private:
    // alias to discoverable data member for ease of notation and backwards compatibility
    const double& _c{_scaleLength};

    // data members initialized during setup
    double _rho0{0.};
};

////////////////////////////////////////////////////////////////////

#endif
