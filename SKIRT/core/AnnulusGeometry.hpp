/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef ANNULUSGEOMETRY_HPP
#define ANNULUSGEOMETRY_HPP

#include "SepAxGeometry.hpp"

////////////////////////////////////////////////////////////////////

/** The AnnulusGeometry class is a subclass of the SepAxGeometry class, and describes a uniform
    cylindrical geometry with a cylindrical hole. This geometry can model geometrically thin
    accretion disks with an inner truncation at \f$R_{\text{min}}\f$, or could be used to build
    multi-zone disks by combining many annuli.

    In formula, this is easily expressed in cylindrical coordinates as \f[ \rho(R, z) = \rho_0\
    \quad\text{for } R_{\text{min}}<R<R_{\text{max}} \text{ and } -\frac{h}{2}<z<\frac{h}{2}. \f]

    The model has three free parameters: the cylinder height \f$h\f$, the radial outer radius
    \f$R_{\text{max}}\f$ and the inner truncation radius \f$R_{\text{min}}\f$. */
class AnnulusGeometry : public SepAxGeometry
{
    ITEM_CONCRETE(AnnulusGeometry, SepAxGeometry, "an annulus geometry")

        PROPERTY_DOUBLE(minRadius, "the inner radius of the annulus")
        ATTRIBUTE_QUANTITY(minRadius, "length")
        ATTRIBUTE_MIN_VALUE(minRadius, "[0")
        ATTRIBUTE_DEFAULT_VALUE(minRadius, "0")

        PROPERTY_DOUBLE(maxRadius, "the outer radius of the annulus")
        ATTRIBUTE_QUANTITY(maxRadius, "length")
        ATTRIBUTE_MIN_VALUE(maxRadius, "]0")

        PROPERTY_DOUBLE(height, "the annulus height")
        ATTRIBUTE_QUANTITY(height, "length")
        ATTRIBUTE_MIN_VALUE(height, "]0")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function verifies the validity of the parameters and calculates the uniform density
        \f$\rho_0\f$. This density is set by the normalisation condition that the total mass equals
        one, i.e. \f[ \rho_0 = \frac{1}{V} = \frac{1}{(R_{\text{max}}^2 - R_{\text{min}}^2) \pi h}.
        \f] */
    void setupSelfBefore() override;

    //======================== Other Functions =======================

    /** This function returns the density \f$\rho(R,z)\f$ at the cylindrical radius \f$R\f$ and
        height \f$z\f$. It just implements the geometry definition. */
    double density(double R, double z) const override;

    /** This function returns the cylindrical radius \f$R\f$ of a random position drawn from the
        geometry, by picking a uniform deviate \f${\cal{X}}\f$ and setting \f$R =
        \sqrt{R_\text{min}^2 + {\cal{X}}_1\,(R_\text{max}^2-R_\text{min}^2)} \f$. */
    double randomCylRadius() const override;

    /** This function returns the height \f$z\f$ of a random position drawn from the geometry, by
        picking a uniform deviate \f${\cal{X}}\f$ and setting \f$z = -\frac{h}{2} + {\cal{X}} h\f$.
        */
    double randomZ() const override;

    /** This function returns the surface density along a line in the equatorial plane starting at
        the centre of the coordinate system, i.e. \f[ \Sigma_R = \int_0^\infty \rho(R,0)\,
        {\text{d}}R. \f] For the annulus geometry, we find \f[ \Sigma_R = \rho_0 (R_{\text{max}} -
        R_{\text{min}}). \f] */
    double SigmaR() const override;

    /** This function returns the surface density along the Z-axis, i.e. the integration of the
        density along the entire Z-axis, \f[ \Sigma_Z = \int_{-\infty}^\infty \rho(0,z)\,
        {\text{d}}z.\f] For the annulus geometry, we find \f[ \Sigma_Z = \rho_0 h. \f] If there is
        an inner hole (\f$R_{\text{min}}>0\f$), obviously \f$\Sigma_Z=0\f$. */
    double SigmaZ() const override;

    //======================== Data Members ========================

private:
    // aliases to discoverable data members for ease of notation and backwards compatibility
    const double& _Rmin{_minRadius};
    const double& _Rmax{_maxRadius};
    const double& _h{_height};

    // data members initialized during setup
    double _rho0{0.};
};

////////////////////////////////////////////////////////////////////

#endif
