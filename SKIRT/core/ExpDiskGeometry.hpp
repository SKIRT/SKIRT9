/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef EXPDISKGEOMETRY_HPP
#define EXPDISKGEOMETRY_HPP

#include "SepAxGeometry.hpp"

////////////////////////////////////////////////////////////////////

/** The ExpDiskGeometry class is a subclass of the SepAxGeometry class, and describes axisymmetric
    geometries characterized by a double-exponential profile, in which the density decreases
    exponentially in the radial and the vertical directions; see van der Kruit (1986, A&A, 157,
    230–244). A truncation can be applied in both the radial and vertical direction, and an inner
    cylindrical hole can be included. In formula form \f[ \rho(R,z) = \rho_0\,
    {\text{e}}^{-\frac{R}{h_R}-\frac{|z|}{h_z}}, \f] for \f$R_{\text{min}} \leq R \leq
    R_{\text{max}}\f$ and \f$|z|\leq z_{\text{max}}\f$. The model contains five free parameters:
    the scale length \f$h_R\f$, the vertical scale height \f$h_z\f$, the radial truncation
    radius \f$R_{\text{max}}\f$, the vertical truncation radius \f$z_{\text{max}}\f$, and the inner
    radius \f$R_{\text{min}}\f$. */
class ExpDiskGeometry : public SepAxGeometry
{
    ITEM_CONCRETE(ExpDiskGeometry, SepAxGeometry, "an exponential disk geometry")

        PROPERTY_DOUBLE(scaleLength, "the scale length")
        ATTRIBUTE_QUANTITY(scaleLength, "length")
        ATTRIBUTE_MIN_VALUE(scaleLength, "]0")

        PROPERTY_DOUBLE(scaleHeight, "the scale height")
        ATTRIBUTE_QUANTITY(scaleHeight, "length")
        ATTRIBUTE_MIN_VALUE(scaleHeight, "]0")

        PROPERTY_DOUBLE(minRadius, "the radius of the central cavity")
        ATTRIBUTE_QUANTITY(minRadius, "length")
        ATTRIBUTE_MIN_VALUE(minRadius, "[0")
        ATTRIBUTE_DEFAULT_VALUE(minRadius, "0")

        PROPERTY_DOUBLE(maxRadius, "the truncation radius (zero means no truncation)")
        ATTRIBUTE_QUANTITY(maxRadius, "length")
        ATTRIBUTE_MIN_VALUE(maxRadius, "[0")
        ATTRIBUTE_DEFAULT_VALUE(maxRadius, "0")
        ATTRIBUTE_DISPLAYED_IF(maxRadius, "Level2")

        PROPERTY_DOUBLE(maxZ, "the truncation height (zero means no truncation)")
        ATTRIBUTE_QUANTITY(maxZ, "length")
        ATTRIBUTE_MIN_VALUE(maxZ, "[0")
        ATTRIBUTE_DEFAULT_VALUE(maxZ, "0")
        ATTRIBUTE_DISPLAYED_IF(maxZ, "Level2")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function verifies the validity of the parameters. The central density \f$\rho_0\f$ is
        set by the normalization condition that the total mass equals one. One finds after some
        elementary calculus \f[ \frac{1}{\rho_0} = 4\pi\, h_R^2\, h_z \left( 1 -
        {\text{e}}^{-z_{\text{max}}/h_z} \right) \left[ \left( 1+\frac{R_{\text{min}}}{h_R} \right)
        {\text{e}}^{-R_{\text{min}}/h_R}- \left( 1+\frac{R_{\text{max}}}{h_R} \right)
        {\text{e}}^{-R_{\text{max}}/h_R} \right] . \f] In case there is no truncation in either
        radial or vertical directions and no inner hole, this reduces to \f[ \rho_0 = \frac{1}{
        4\pi\, h_R^2\, h_z }. \f] */
    void setupSelfBefore() override;

    //======================== Other Functions =======================

    /** This function returns the density \f$\rho(R,z)\f$ at the cylindrical radius \f$R\f$ and
        height \f$z\f$. It just implements the analytical formula. */
    double density(double R, double z) const override;

    /** This function returns the cylindrical radius \f$R\f$ of a random position drawn from the
        geometry, by picking a uniform deviate \f${\cal{X}}\f$ and solving the equation \f[
        {\cal{X}} = 2\pi \int_0^R \rho_R(R')\, R'\, {\text{d}}R' \f] for \f$R\f$. Substituting the
        exponential radial profile (without truncation) into this equation, we obtain \f[ {\cal{X}}
        = 1 - \left( 1+\frac{R}{h_R} \right) \exp \left( -\frac{R}{h_R} \right). \f] This equation
        can be solved by means of the Lambert function of order \f$-1\f$, yielding \f[ R = h_R
        \left[ -1-W_{-1} \left( \frac{ {\cal{X}}-1}{\text{e}} \right) \right]. \f] The Lambert
        function \f$W_{-1}(z)\f$ is implemented in the function SpecialFunctions::LambertW1. The
        truncation and the inner hole are taken into account by rejecting values larger than
        \f$R_{\text{max}}\f$ or smaller than \f$R_{\text{min}}\f$. */
    double randomCylRadius() const override;

    /** This function returns the height \f$z\f$ of a random position drawn from the geometry, by
        picking a uniform deviate \f${\cal{X}}\f$ and solving the equation \f[ {\cal{X}} =
        \int_{-\infty}^z \rho_z(z')\, {\text{d}}z' \f] for \f$z\f$. For the exponential disk
        geometry, this integration is simple, and the inversion results in \f[ z = \begin{cases} \;
        h_z\,\ln(2{\cal{X}}) & \text{if $0<{\cal{X}}<\tfrac{1}{2}$,} \\ \;-h_z\,\ln[2(1-{\cal{X}})]
        & \text{if $\tfrac{1}{2}<{\cal{X}}<1$.} \end{cases} \f] The truncation is taken into
        account by rejecting values \f$|z|\f$ larger than \f$z_{\text{max}}\f$. */
    double randomZ() const override;

    /** This function returns the surface density along a line in the equatorial plane starting at
        the centre of the coordinate system, i.e. \f[ \Sigma_R = \int_0\infty \rho(R,0)\,
        {\text{d}}R. \f] For the exponential disc geometry we find \f[ \Sigma_R = \rho_0 h_R \left(
        {\text{e}}^{-R_{\text{min}}/h_R} - {\text{e}}^{-R_{\text{max}}/h_R} \right), \f] which
        reduces to \f$ \Sigma_R = \rho_0 h_R \f$ if there is no radial truncation and no inner
        hole. */
    double SigmaR() const override;

    /** This function returns the surface density along the Z-axis, i.e. the integration of the
        density along the entire Z-axis, \f[ \Sigma_Z = \int_{-\infty}^\infty \rho(0,0,z)\,
        {\text{d}}z.\f] For the exponential disc geometry we find \f[ \Sigma_Z = 2\,\rho_0 h_Z
        \left( 1 - {\text{e}}^{-z_{\text{max}}/h_z} \right), \f] which reduces to \f$ \Sigma_Z =
        2\,\rho_0 h_z \f$ if there is no vertical truncation. If there is an inner hole, obviously
        \f$\Sigma_Z=0\f$. */
    double SigmaZ() const override;

    //======================== Data Members ========================

private:
    // aliases to discoverable data members for ease of notation and backwards compatibility
    const double& _hR{_scaleLength};
    const double& _hz{_scaleHeight};
    const double& _Rmin{_minRadius};
    const double& _Rmax{_maxRadius};
    const double& _zmax{_maxZ};

    // data members initialized during setup
    double _rho0{0.};
};

////////////////////////////////////////////////////////////////////

#endif
