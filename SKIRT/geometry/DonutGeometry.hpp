/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef DONUTGEOMETRY_HPP
#define DONUTGEOMETRY_HPP

#include "Array.hpp"
#include "AxGeometry.hpp"

////////////////////////////////////////////////////////////////////

/** The DonutGeometry class is a subclass of the AxGeometry class and describes the geometry of an
    axisymmetric ring torus (i.e. a doughnut geometry) of uniform density. This ring torus geometry
    is often used to model the dusty gas in the equatorial plane of active galactic nuclei (AGN),
    to explain the observed variety of AGN types as a distribution of observed system orientations
    relative to this obscuring torus of gas and dust; see Antonucci (1993, ARA&A, 31, 473), Urry &
    Padovani (1995, PASP, 107, 803) and Netzer (2015, ARA&A, 53, 365).

    The ring torus geometry is obtained by revolving a circle of radius \f$r_{\text{small}}\f$
    about an axis that is co-planar with this circle, at a radius \f$r_{\text{large}}\f$. Inside
    this surface of revolution the density has a constant value, and outside this surface the
    density is zero. In formula, this is easily expressed in cylindrical coordinates as \f[ \rho(R,
    z) = A\ \quad\text{for } \sqrt{(R-r_{\text{large}})^2 + z^2}<r_{\text{small}}. \f] There are
    two free parameters describing this torus geometry: the torus radius in the azimuthal plane
    \f$r_{\text{small}}\f$ and the torus radius in the equatorial plane \f$r_{\text{large}}\f$. */
class DonutGeometry : public AxGeometry
{
    ITEM_CONCRETE(DonutGeometry, AxGeometry, "a donut torus geometry")

        PROPERTY_DOUBLE(largeRadius, "the radius of the torus circle in the equatorial plane")
        ATTRIBUTE_QUANTITY(largeRadius, "length")
        ATTRIBUTE_MIN_VALUE(largeRadius, "]0")

        PROPERTY_DOUBLE(smallRadius, "the radius of the torus circle in the azimuthal plane")
        ATTRIBUTE_QUANTITY(smallRadius, "length")
        ATTRIBUTE_MIN_VALUE(smallRadius, "]0")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function calculates some frequently used values. The normalization parameter \f$A\f$
        is set by the normalization condition that total mass equals one, i.e. \f[ A = \frac{1}{V}
        = \frac{1}{2 \pi^2 r_{\text{large}} r_{\text{small}}^2}. \f] It also tabulates, using
        trapezoidal integration at a large number of radii spanning \f$[r_{\text{large}} -
        r_{\text{small}}, r_{\text{large}} + r_{\text{small}}]\f$, the cumulative distribution of
        the marginal radial probability density \f[ p(R)\,{\text{d}}R \propto R\,
        \sqrt{r_{\text{small}}^2 - (R-r_{\text{large}})^2}\,{\text{d}}R, \f] so that
        generatePosition() can generate a random radius through numerical inversion. */
    void setupSelfBefore() override;

    //======================== Other Functions =======================

public:
    /** This function returns the density \f$\rho(R,z)\f$ at the cylindrical radius \f$R\f$ and
        height \f$z\f$. It just implements the geometry definition. */
    double density(double R, double z) const override;

    /** This function generates a random position from the uniform ring torus geometry. Because
        the cylindrical volume element \f$R\,{\text{d}}R\,{\text{d}}\varphi\,{\text{d}}z\f$
        depends on \f$R\f$, sampling a position uniformly within the meridional cross-section
        circle and revolving it around the z-axis would not produce a position that is uniform by
        volume within the torus. Instead, the function first draws a random cylindrical radius
        \f$R\f$ from the marginal probability density \f$p(R)\,{\text{d}}R \propto R\,
        \sqrt{r_{\text{small}}^2-(R-r_{\text{large}})^2}\,{\text{d}}R\f$ tabulated by
        setupSelfBefore(), using linear interpolation of its cumulative distribution (the
        cumulative distribution of this \f$p(R)\f$ does not invert analytically, because its
        antiderivative mixes an arcsine term with algebraic terms in \f$R\f$). Because the density
        is uniform inside the torus, the height \f$z\f$ then has a uniform conditional
        distribution given this \f$R\f$, over the range \f$|z| \le
        \sqrt{r_{\text{small}}^2-(R-r_{\text{large}})^2}\f$ set by the torus tube at that radius,
        so a random \f$z\f$ is drawn directly by choosing a random deviate \f${\cal{X}}_1\f$ and
        setting \f$z = \sqrt{r_{\text{small}}^2-(R-r_{\text{large}})^2}\,(2{\cal{X}}_1-1)\f$.
        Finally, a random azimuth angle \f$\phi\f$ is found by choosing a random deviate
        \f${\cal{X}}_2\f$ and setting \f$\phi = 2\pi {\cal{X}}_2\f$. */
    Position generatePosition() const override;

    /** This function returns the radial surface density, i.e. the integration of the density along
        a line in the equatorial plane starting at the centre of the coordinate system, \f[
        \Sigma_R = \int_0^\infty \rho(R,0)\,{\text{d}}R. \f] For this geometry, \f[ \Sigma_R = A\,2
        r_{\text{small}}. \f] */
    double SigmaR() const override;

    /** This function returns the Z-axis surface density, i.e. the integration of the density along
        the entire Z-axis, \f[ \Sigma_Z = \int_{-\infty}^\infty \rho(0,z)\, {\text{d}}z. \f] For
        this geometry, the integral is simply zero. */
    double SigmaZ() const override;

    //======================== Data Members ========================

private:
    // aliases to discoverable data members for ease of notation and backwards compatibility
    const double& _R0{_largeRadius};
    const double& _r0{_smallRadius};

    // data members initialized during setup
    double _A{0.};
    double _r02{0.};
    Array _Rv;
    Array _Xv;
};

////////////////////////////////////////////////////////////////////

#endif
