/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef SHELLGEOMETRY_HPP
#define SHELLGEOMETRY_HPP

#include "SpheGeometry.hpp"

////////////////////////////////////////////////////////////////////

/** The ShellGeometry class is a subclass of the SpheGeometry class and describes the geometry of a
    spherical shell, where the density behaves as a power law between an inner and an outer radius,
    \f[ \rho(r) = A\,r^{-p} \qquad\qquad r_{\text{min}} < r < r_{\text{max}}. \f] with \f$A\f$ a
    normalization constant. The range of \f$p\f$ is limited to \f$p\geq0\f$, and obviously the
    condition \f$r_{\text{min}} < r_{\text{max}}\f$ should be satisfied. This geometry is
    characterized by three free parameters: the inner radius \f$r_{\text{min}}\f$, the outer radius
    \f$r_{\text{max}}\f$ and the power law exponent \f$p\f$. */
class ShellGeometry : public SpheGeometry
{
    ITEM_CONCRETE(ShellGeometry, SpheGeometry, "a shell geometry")

        PROPERTY_DOUBLE(minRadius, "the inner radius of the shell")
        ATTRIBUTE_QUANTITY(minRadius, "length")
        ATTRIBUTE_MIN_VALUE(minRadius, "]0")

        PROPERTY_DOUBLE(maxRadius, "the outer radius of the shell")
        ATTRIBUTE_QUANTITY(maxRadius, "length")
        ATTRIBUTE_MIN_VALUE(maxRadius, "]0")

        PROPERTY_DOUBLE(exponent, "the power law exponent")
        ATTRIBUTE_MIN_VALUE(exponent, "[0")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function verifies the validity of the attributes. The normalization parameter \f$A\f$
        is set by the normalization condition that total mass equals one, i.e. \f[ 1 = 4\pi A
        \int_{r_{\text{min}}}^{r_{\text{max}}} r^{2-p}\, {\text{d}}r. \f] This results in \f[ A =
        \frac{1}{4\pi}\, \frac{1}{ {\text{gln}}_{p-2}\, r_{\text{max}} - {\text{gln}}_{p-2}\,
        r_{\text{min}} }, \f] with \f${\text{gln}}_p\,x\f$ the generalized logarithm defined in
        SpecialFunctions::gln. */
    void setupSelfBefore() override;

    //======================== Other Functions =======================

public:
    /** This function returns the density \f$\rho(r)\f$ at the radius \f$r\f$. It just implements
        the analytical formula. */
    double density(double r) const override;

    /** This function returns the radius of a random position drawn from the shell density
        distribution. This is accomplished by generating a uniform deviate \f${\cal{X}}\f$, and
        solving the equation \f[ {\cal{X}} = M(r) = 4\pi \int_0^r \rho(r')\, r'{}^2\, {\text{d}}r'
        \f] for \f$r\f$. For the shell geometry, with \f$p\ne3\f$, we find \f[ {\cal{X}} =
        \frac{r^{3-p}-r_{\text{min}}^{3-p}} {r_{\text{max}}^{3-p}-r_{\text{min}}^{3-p}}. \f]
        Inverting this results in \f[ r = \left[ (1-{\cal{X}})\,r_{\text{min}}^{3-p} +
        {\cal{X}}\,r_{\text{max}}^{3-p} \right]^{\frac{1}{3-p}}. \f] For \f$p=3\f$ this expression
        does not hold, and for \f$p\approx3\f$ it breaks down numerically. So for \f$p\approx3\f$
        we can write the general expression \f[ r = {\text{gexp}}_{p-2} \Big[ {\text{gln}}_{p-2}\,
        r_{\text{min}} + {\cal{X}}\,( {\text{gln}}_{p-2}\, r_{\text{max}} - {\text{gln}}_{p-2}\,
        r_{\text{min}} ) \Bigr]. \f] In this expression, \f${\text{gln}}_p\,x\f$ and
        \f${\text{gexp}}_p\,x\f$ are the generalized logarithm and exponential functions defined in
        SpecialFunctions::gln and SpecialFunctions::gexp respectively. */
    double randomRadius() const override;

    /** This function returns the radial surface density, i.e. the integration of the density along
        a line starting at the centre of the coordinate system, \f[ \Sigma_r = \int_0^\infty
        \rho(r)\,{\text{d}}r. \f] For the shell geometry, one obtains \f[ \Sigma_r = A\, (
        {\text{gln}}_p\, r_{\text{max}} - {\text{gln}}_p\, r_{\text{min}} ) \f] with
        \f${\text{gln}}_p\,x\f$ the generalized logarithm defined in SpecialFunctions::gln. */
    double Sigmar() const override;

    //======================== Data Members ========================

private:
    // aliases to discoverable data members for ease of notation and backwards compatibility
    const double& _rmin{_minRadius};
    const double& _rmax{_maxRadius};
    const double& _p{_exponent};

    // data members initialized during setup
    double _smin{0.};
    double _sdiff{0.};
    double _tmin{0.};
    double _tmax{0.};
    double _A{0.};
};

////////////////////////////////////////////////////////////////////

#endif
