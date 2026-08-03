/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef RINGGEOMETRY_HPP
#define RINGGEOMETRY_HPP

#include "Array.hpp"
#include "SepAxGeometry.hpp"

//////////////////////////////////////////////////////////////////////

/** The RingGeometry class is a subclass of the SepAxGeometry class and describes the geometry
    of a ring as one would consider in dust-lane early-type galaxies. The geometry is
    characterized by a gaussian profile in the radial direction and an exponential fall-off in the
    vertical direction, resulting in the axisymmetric density profile \f[ \rho(R,z) = A\,
    \exp\left[ -\frac{(R-R_0)^2}{2w^2} \right] \exp\left( -\frac{|z|}{h_z} \right). \f] There are
    three free parameters in this class: the radius \f$R_0\f$ of the ring, the radial dispersion
    \f$w\f$ of the ring and the vertical scale height \f$h_z\f$. */
class RingGeometry : public SepAxGeometry
{
    ITEM_CONCRETE(RingGeometry, SepAxGeometry, "a ring geometry")

        PROPERTY_DOUBLE(ringRadius, "the radius of the ring")
        ATTRIBUTE_QUANTITY(ringRadius, "length")
        ATTRIBUTE_MIN_VALUE(ringRadius, "]0")

        PROPERTY_DOUBLE(width, "the radial width of the ring")
        ATTRIBUTE_QUANTITY(width, "length")
        ATTRIBUTE_MIN_VALUE(width, "]0")

        PROPERTY_DOUBLE(height, "the scale height of the ring")
        ATTRIBUTE_QUANTITY(height, "length")
        ATTRIBUTE_MIN_VALUE(height, "]0")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function verifies the validity of the parameters of the geometry. The density
        normalization parameter \f$A\f$ is set by the normalization condition that the total
        mass equals one. We find \f[ A = \frac{1}{4\pi\,w^2\,h_z
        \left[ {\text{e}}^{-t^2} + \sqrt{\pi}\,t\,(1 + {\text{erf}}\,t) \right]}, \f] with
        \f$t = R_0/\sqrt{2}w\f$. The routine also creates a vector with the cumulative distribution
        \f$X(r)\f$ corresponding to the radial part of the density, i.e.
        mass \f[ X(R) = \int_0^R \rho_R(R')\, R'\, {\text{d}}R' \f] at a large number
        of radii. For the ring geometry, \f[ \rho_R(R') = 4\pi\,A\,h_z \exp\left[
        -\frac{(R-R_0)^2}{2w^2} \right], \f] and thus \f[ X(R) =
        4\pi\,A\,w^2\,h_z \left[ \left({\text{e}}^{-t^2}-{\text{e}}^{-u^2}\right)
        + \sqrt{\pi}\,t\,({\text{erf}}\,t - {\text{erf}}\,u) \right], \f] with \f$t\f$ as defined before
        and \f$u = (R_0-R)/\sqrt{2}w\f$. */
    void setupSelfBefore() override;

    //======================== Other Functions =======================

public:
    /** This function returns the density \f$\rho(R,z)\f$ at the cylindrical radius
        \f$R\f$ and height \f$z\f$. It just implements the analytical formula. */
    double density(double R, double z) const override;

    /** This function returns the cylindrical radius \f$R\f$ of a random position drawn from the
        geometry, by picking a uniform deviate \f${\cal{X}}\f$
        and solving the equation \f[ {\cal{X}} = 2\pi \int_0^R \rho_R(R')\, R'\, {\text{d}}R' \f]
        for \f$R\f$. We just use the vector of cumulative masses stored internally. */
    double randomCylRadius() const override;

    /** This function returns the height \f$z\f$ of a random position drawn from the
        geometry, by picking a uniform deviate \f${\cal{X}}\f$
        and solving the equation \f[ {\cal{X}} = \int_{-\infty}^z \rho_z(z')\, {\text{d}}z' \f]
        for \f$z\f$. For the ring geometry with its exponential vertical profile, this
        integration is simple, and the inversion results in \f[ z = \begin{cases} \;
        h_z\,\ln(2{\cal{X}}) & \text{if $0<{\cal{X}}<\tfrac{1}{2}$,} \\ \;-h_z\,\ln[2(1-{\cal{X}})]
        & \text{if $\tfrac{1}{2}<{\cal{X}}<1$.} \end{cases} \f] */
    double randomZ() const override;

    /** This function returns the surface density along a line in the equatorial plane
        starting at the centre of the coordinate system, i.e. \f[ \Sigma_R =
        \int_0^infty \rho(R,0)\,{\text{d}}R. \f] For the
        ring geometry we obtain \f[ \Sigma_R = \sqrt{\frac{\pi}{2}}\,A\,w \left[ 1 +
        {\text{erf}} \left( \frac{R_0}{\sqrt{2} w} \right) \right]. \f] */
    double SigmaR() const override;

    /** This function returns the surface density along the Z-axis, i.e.
        the integration of the density along the entire Z-axis, \f[
        \Sigma_Z = \int_{-\infty}^\infty \rho(0,0,z)\, {\text{d}}z.\f] For the ring geometry we
        find \f[ \Sigma_{\text{f}} = 2A\,h_z\, \exp\left(-\frac{R_0^2}{2w^2}\right). \f] */
    double SigmaZ() const override;

    //======================== Data Members ========================

private:
    // aliases to discoverable data members for ease of notation and backwards compatibility
    const double& _R0{_ringRadius};
    const double& _w{_width};
    const double& _hz{_height};

    // data members initialized during setup
    double _A{0.};
    Array _Rv;
    Array _Xv;
};

//////////////////////////////////////////////////////////////////////

#endif
