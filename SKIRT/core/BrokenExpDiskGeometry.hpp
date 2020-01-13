/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef BROKENEXPDISKGEOMETRY_HPP
#define BROKENEXPDISKGEOMETRY_HPP

#include "Array.hpp"
#include "SepAxGeometry.hpp"

////////////////////////////////////////////////////////////////////

/** The BrokenExpDiskGeometry class is a subclass of the SepAxGeometry class. It describes a
    particular class of models for discs with a so-called break in the radial profile. More
    specifically, the radial density profile is a broken exponential profile: an inner exponential
    profile with scale length \f$h_{R,{\text{inn}}}\f$ and an outer exponential profile with scale
    length \f$h_{R,{\text{out}}}\f$. There is no limitation in the nature of the break: the inner
    profile can be shallower than the outer one (\f$h_{R,{\text{out}}} < h_{R,{\text{inn}}}\f$) or
    steeper (\f$h_{R,{\text{out}}} > h_{R,{\text{inn}}}\f$). These two asymptotic profiles are
    joined smoothly at a break radius \f$R_{\text{b}}\f$, with a dimensionless parameter \f$s\f$
    that sets the smoothness of the transition between these two regimes. For large values of the
    smoothness parameter, \f$s\gg1\f$, the transition is sharp, whereas for small values,
    \f$s\sim1\f$, the transition is very gradual. In the vertical direction, the density decreases
    exponentially. No truncation is applied. In formula form, \f[ \rho(R,z) = \rho_0
    {\text{e}}^{-\frac{R}{h_{\text{inn}}}-\frac{|z|}{h_z}} \left( 1 +
    {\text{e}}^{\frac{s\,(R-R_{\text{b}})}{h_{\text{out}}}} \right)^{\frac{1}{s}
    \left(\frac{h_{\text{out}}}{h_{\text{inn}}} - 1\right)}. \f] This geometry is strongly inspired
    by Erwin et al. (2008, AJ, 135, 20) and Erwin (2015, ApJ, 799, 226), but has one significant
    difference: in our formulation, the sharpness \f$s\f$ is a dimensionless number, whereas in
    Erwin et al. (2008) the sharpness parameter \f$\alpha\f$ is a quantity with dimension
    length\f$^{-1}\f$. More specifically, \f[ \alpha = \frac{s}{h_{\text{out}}} \f] The model
    contains five free parameters: the scale length of inner disc \f$h_{\text{inn}}\f$, the scale
    length of outer disc \f$h_{\text{out}}\f$, the scale height \f$h_z\f$, the break radius
    \f$R_{\text{b}}\f$, and the sharpness of the break \f$s\f$. The final parameter that appears in
    the formula above is the central density \f$\rho_0\f$; it is not a free parameter, but such
    that the total mass of the geometry is normalized to one. */
class BrokenExpDiskGeometry : public SepAxGeometry
{
    ITEM_CONCRETE(BrokenExpDiskGeometry, SepAxGeometry, "a broken exponential disk geometry")
        ATTRIBUTE_TYPE_DISPLAYED_IF(BrokenExpDiskGeometry, "Level2")

        PROPERTY_DOUBLE(scaleLengthInner, "the scale length of the inner disk")
        ATTRIBUTE_QUANTITY(scaleLengthInner, "length")
        ATTRIBUTE_MIN_VALUE(scaleLengthInner, "]0")

        PROPERTY_DOUBLE(scaleLengthOuter, "the scale length of the outer disk")
        ATTRIBUTE_QUANTITY(scaleLengthOuter, "length")
        ATTRIBUTE_MIN_VALUE(scaleLengthOuter, "]0")

        PROPERTY_DOUBLE(scaleHeight, "the scale height")
        ATTRIBUTE_QUANTITY(scaleHeight, "length")
        ATTRIBUTE_MIN_VALUE(scaleHeight, "]0")

        PROPERTY_DOUBLE(breakRadius, "the break radius")
        ATTRIBUTE_QUANTITY(breakRadius, "length")
        ATTRIBUTE_MIN_VALUE(breakRadius, "]0")

        PROPERTY_DOUBLE(sharpness, "the sharpness of the break")
        ATTRIBUTE_MIN_VALUE(sharpness, "]0")
        ATTRIBUTE_DEFAULT_VALUE(sharpness, "3")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function calculates some frequently used values. The central density \f$\rho_0\f$ is
        set by the normalization condition that the total mass equals one. One finds after some
        elementary calculus \f[ \rho_0 = \frac{1}{4\pi\,h_z\, I_R} \f] where \f$ I_R\f$ is the
        integral \f[ I_R = \int_0^\infty \rho_R(R)\,R\, {\text{d}}R \f] with \f$\rho_R(R)\f$ the
        radial part of the density distribution. This integral is calculated using a simple
        trapezoidal integration. This routine immediately also stores the cumulative radial
        distribution \f[ X(R) = \frac{1}{I_R} \int_0^R \rho_R(R')\, R'\, {\text{d}}R' \f] in an
        internal array, so that this can be used later on to generate random positions extracted
        from this geometry. Finally, also the radial surface density is calculated using a similar
        numerical integration. */
    void setupSelfBefore() override;

    //======================== Other Functions =======================

    /** This function returns the density \f$\rho(R,z)\f$ at the cylindrical radius \f$R\f$ and
        height \f$z\f$. It just implements the analytical formula. */
    double density(double R, double z) const override;

    /** This function returns the cylindrical radius \f$R\f$ of a random position drawn from the
        geometry, by picking a uniform deviate \f${\cal{X}}\f$. We just use the vector of
        cumulative masses stored internally. */
    double randomCylRadius() const override;

    /** This function returns the height \f$z\f$ of a random position drawn from the geometry, by
        picking a uniform deviate \f${\cal{X}}\f$ and solving the equation \f[ {\cal{X}} =
        \int_{-\infty}^z \rho_z(z')\, {\text{d}}z' \f] for \f$z\f$. For the exponential disk
        geometry, this integration is simple, and the inversion results in \f[ z = \begin{cases} \;
        h_z\,\ln(2{\cal{X}}) & \text{if $0<{\cal{X}}<\tfrac{1}{2}$,} \\ \;-h_z\,\ln[2(1-{\cal{X}})]
        & \text{if $\tfrac{1}{2}<{\cal{X}}<1$.} \end{cases} \f] */
    double randomZ() const override;

    /** This function returns the surface density along a line in the equatorial plane starting at
        the centre of the coordinate system, i.e. \f[ \Sigma_R = \int_0^\infty \rho(R,z)\,
        {\text{d}}R. \f] This value is calculated numerically during setup and stored as a data
        member. */
    double SigmaR() const override;

    /** This function returns the surface density along the Z-axis, i.e. the integration of the
        density along the entire Z-axis, \f[ \Sigma_Z = \int_{-\infty}^\infty \rho(0,z)\,
        {\text{d}} z \f] For the BrokenExpDiskGeometry model, one obtains \f[ \Sigma_Z = 2\,
        \rho_0\, z_0 \left(1+{\text{e}}^{-s\,R_{break}}\right)^{\frac{1}{s}
        \left(\frac{h_{\text{out}}}{h_{\text{inn}}} - 1\right)}\f] */
    double SigmaZ() const override;

private:
    /** This private member function returns the radial part of the density distribution, i.e. \f[
        \rho_R(R) = {\text{e}}^{-\frac{R}{h_{\text{inn}}}} \left( 1 +
        {\text{e}}^{\frac{s\,(R-R_{\text{b}})}{h_{\text{out}}}} \right)^{\frac{1}{s}
        \left(\frac{h_{\text{out}}}{h_{\text{inn}}} - 1\right)} \f] It is just an auxiliary
        function. */
    double radialDensity(double R) const;

    //======================== Data Members ========================

private:
    // aliases to discoverable data members for ease of notation and backwards compatibility
    const double& _hinn{_scaleLengthInner};
    const double& _hout{_scaleLengthOuter};
    const double& _hz{_scaleHeight};
    const double& _Rb{_breakRadius};
    const double& _s{_sharpness};

    // data members initialized during setup
    double _beta{0.};
    double _rho0{0.};
    double _SigmaR{0.};
    Array _Rv;
    Array _Xv;
};

////////////////////////////////////////////////////////////////////

#endif
