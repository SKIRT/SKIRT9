/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef TTAURIDISKGEOMETRY_HPP
#define TTAURIDISKGEOMETRY_HPP

#include "AxGeometry.hpp"

////////////////////////////////////////////////////////////////////

/** The TTauriDiskGeometry class describes the geometry of a typical passive disk around a T Tauri
    star. The disks are axisymmetric with a central cavity and are characterized by the density
    distribution

    \f[ \rho(R,z) = \rho_0 \left(\frac{R}{R_d}\right)^{-\alpha} \exp\left\{ -b \left[ \frac{z/z_d}
    {(R/R_d)^{9/8}} \right]^2 \right\} \qquad\qquad R_{\text{inn}} < R < R_{\text{out}}. \f]

    There are six parameters for this geometry: the exponent indices \f$\alpha>0\f$ and \f$b>0\f$,
    the inner and outer radii \f$R_{\text{inn}}>0\f$ and \f$R_{\text{out}}>R_{\text{inn}}\f$, and
    the scale length \f$R_d>0\f$ and scale height \f$z_d>0\f$.

    Special cases of this geometry are used in the radiative transfer benchmark problems described
    by Pascucci et al. (2004, A&A, 417, 793) and Pinte et al. (2009, A&A, 498, 967).

    Note: for the implementation of this geometry we use the generalized logarithm defined by
    SpecialFunctions::gln() and its inverse SpecialFunctions::gexp(). Specifically, we often use
    the difference function SpecialFunctions::gln2() defined as \f[ {\text{gln2}}(p,x_1,x_2) =
    {\text{gln}}(p,x_1) - {\text{gln}}(p,x_2) = \int_{x_2}^{x_1} t^{-p}\,{\text{d}}t. \f] */
class TTauriDiskGeometry : public AxGeometry
{
    ITEM_CONCRETE(TTauriDiskGeometry, AxGeometry, "a T Tauri disk geometry")
        ATTRIBUTE_TYPE_DISPLAYED_IF(TTauriDiskGeometry, "Level2")

        PROPERTY_DOUBLE(scaleLength, "the scale length")
        ATTRIBUTE_QUANTITY(scaleLength, "length")
        ATTRIBUTE_MIN_VALUE(scaleLength, "]0")

        PROPERTY_DOUBLE(scaleHeight, "the scale height")
        ATTRIBUTE_QUANTITY(scaleHeight, "length")
        ATTRIBUTE_MIN_VALUE(scaleHeight, "]0")

        PROPERTY_DOUBLE(minRadius, "the inner radius of the disk")
        ATTRIBUTE_QUANTITY(minRadius, "length")
        ATTRIBUTE_MIN_VALUE(minRadius, "]0")

        PROPERTY_DOUBLE(maxRadius, "the outer radius of the disk")
        ATTRIBUTE_QUANTITY(maxRadius, "length")
        ATTRIBUTE_MIN_VALUE(maxRadius, "]0")

        PROPERTY_DOUBLE(radialIndex, "the radial exponent index α")
        ATTRIBUTE_MIN_VALUE(radialIndex, "]0")
        ATTRIBUTE_DEFAULT_VALUE(radialIndex, "2.5")

        PROPERTY_DOUBLE(verticalIndex, "the vertical exponent index b")
        ATTRIBUTE_MIN_VALUE(verticalIndex, "]0")
        ATTRIBUTE_DEFAULT_VALUE(verticalIndex, "0.5")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function verifies that the parameters are valid and it calculates the normalization
        constant \f$\rho_0\f$ determined by requiring that total mass equals one: \f[ 1 = 2\pi\,
        \int_{R_\text{inn}}^{R_\text{out}} R\,{\text{d}}R \, \int_{-\infty}^\infty
        \rho(R,z)\,{\text{d}}z. \f]

        Given that \f[ \int_{-\infty}^{+\infty} \exp(-a^2x^2)\,\text{d}x = \frac{\sqrt{\pi}}{a}, \f]
        and using the definition of gln2() noted in the class header, we find after some algebra
        that \f[ \rho_0 = \left[ 2\, \pi^{3/2}\, b^{-1/2}\, R_d^2\, z_d\,
        \text{gln2}\left(\alpha-\frac{17}{8}, \frac{R_{\text{out}}}{R_d},
        \frac{R_{\text{inn}}}{R_d}\right) \right]^{-1}. \f] */
    void setupSelfBefore() override;

    //======================== Other Functions =======================

public:
    /** This function returns the density \f$\rho(R,z)\f$ at the cylindrical radius
        \f$R\f$ and height \f$z\f$. It just implements the analytical formula. */
    double density(double R, double z) const override;

    /** This function generates a random position from the geometry by drawing a random point from
        the three-dimensional probability density \f$p({\bf{r}})\, {\text{d}}{\bf{r}} =
        \rho({\bf{r}})\, {\text{d}}{\bf{r}}\f$. In the present case, we accomplish this by picking
        a random cylindrical radius \f$R\f$, a random height \f$z\f$, and a random azimuth
        \f$\phi\f$ from the appropriate one-dimensional probability distribution functions.

        We first generate a random radius \f$R\f$ from the marginal distribution \f[
        p(R)\,{\text{d}}R = 2\pi R\, {\text{d}}R\, \int_{-\infty}^\infty \rho(R,z)\, {\text{d}}z.
        \f] The cumulative distribution is \f[ P(R) = 2\pi\, \int_{R_\text{inn}}^R R'\,{\text{d}}R'
        \, \int_{-\infty}^\infty \rho(R',z)\,{\text{d}}z. \f] Using the definitions of gln() and
        gln2() noted in the class header, the cumulative distribution can be written as \f[ P(R) =
        \frac{\text{gln2}\left(\alpha-\frac{17}{8}, \frac{R}{R_d},
        \frac{R_{\text{inn}}}{R_d}\right)} {\text{gln2}\left(\alpha-\frac{17}{8},
        \frac{R_{\text{out}}}{R_d}, \frac{R_{\text{inn}}}{R_d}\right)} =
        \frac{\text{gln}\left(\alpha-\frac{17}{8}, \frac{R}{R_d}, \right) -
        \text{gln}\left(\alpha-\frac{17}{8}, \frac{R_{\text{inn}}}{R_d}\right)}
        {\text{gln2}\left(\alpha-\frac{17}{8}, \frac{R_{\text{out}}}{R_d},
        \frac{R_{\text{inn}}}{R_d}\right)}. \f] A random \f$R\f$ is generated by picking a uniform
        deviate \f${\cal{X}}\f$ and setting \f${\cal{X}}=P(R)\f$. After resolving for \f$R\f$, and
        using the inverse generalized logarithm gexp(), we obtain \f[ R = R_d\,
        \text{gexp}\left[\alpha-\frac{17}{8}, \text{gln}\left(\alpha-\frac{17}{8},
        \frac{R_{\text{inn}}}{R_d}\right) + {\cal{X}}\, \text{gln2}\left(\alpha-\frac{17}{8},
        \frac{R_{\text{out}}}{R_d}, \frac{R_{\text{inn}}}{R_d}\right) \right]. \f]

        We then generate a random height from the conditional distribution function \f[
        p(z)\,{\text{d}}z = \dfrac{\rho(R,z)\,{\text{d}}z}{\int_{-\infty}^\infty \rho(R,z')\,
        {\text{d}}z'}, \f] where \f$R\f$ is the random cylindrical radius generated before. One
        easily finds that this distribution is a Gaussian distribution with mean zero and
        dispersion \f[ \sigma(R) = \frac{1}{\sqrt{2b}}\,z_d\,\left(\frac{R}{R_d}\right)^{9/8}. \f]
        Generating a random \f$z\f$ from this distribution is easy as the Random class contains a
        gaussian random number generating function.

        Finally, we simply generate the azimuth from a uniform distribution between 0 and
        \f$2\pi\f$. */
    Position generatePosition() const override;

    /** This function returns the radial surface density, i.e. the integration of the density along
        a line in the equatorial plane starting at the centre of the coordinate system, \f[
        \Sigma_R = \int_0^\infty \rho(R,0)\,{\text{d}}R. \f] Using the definition of gln2() noted
        in the class header, we easily find for the T Tauri disk geometry that \f[ \Sigma_R =
        \rho_0\, R_d\, \text{gln2}\left( \alpha, \frac{R_{\text{out}}}{R_d},
        \frac{R_{\text{inn}}}{R_d} \right). \f] */
    double SigmaR() const override;

    /** This function returns the vertical surface density, i.e. the integration of the density
        along the entire Z-axis, \f[ \Sigma_Z = \int_{-\infty}^\infty \rho(0,z)\, {\text{d}}z. \f]
        For the T Tauri disk geometry with its central cylindrical cavity, this integral is zero.
        */
    double SigmaZ() const override;

    //======================== Data Members ========================

private:
    // aliases to discoverable data members for ease of notation
    const double& _Rd{_scaleLength};
    const double& _zd{_scaleHeight};
    const double& _Rinn{_minRadius};
    const double& _Rout{_maxRadius};
    const double& _a{_radialIndex};
    const double& _b{_verticalIndex};

    // data members initialized during setup
    double _a178{0.};
    double _glnInn{0.};
    double _glnInnOut{0.};
    double _rho0{0.};
    double _s0{0.};
};

////////////////////////////////////////////////////////////////////

#endif
