/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef HYPERBOLOIDSHELLGEOMETRY_HPP
#define HYPERBOLOIDSHELLGEOMETRY_HPP

#include "AxGeometry.hpp"

////////////////////////////////////////////////////////////////////

/** The HyperboloidShellGeometry class is a subclass of the AxGeometry class and describes the
    geometry defined by two hyperboloid surfaces. It may be used to represent dusty outflows,
    extending vertically and then flaring towards the outside, asymptotically approaching cone
    surface.

    \image html HyperboloidShellGeometry.png

    The current implementation allows only a constant density distribution bounded by two
    hyperboloid surfaces of the form \f[ z=\frac{c}{a}\sqrt{\rho^2-a^2} \f] and two horizontal
    planes paralel to the xy plane, \f$\pm z_{\text{max}}\f$. There are five free parameters
    describing this dust geometry: the radial extent of the outer wall \f$D_\text{in}\f$, the half
    opening angles of the outer and inner walls (\f$\Delta_\text{out}, \Delta_\text{in}\f$), the
    real axes of the outer and inner wall (\f$a_\text{out}, a_\text{in}\f$). The radial extent and
    half opening angles correspond to the lenght of the side and half opening angle of the cone to
    which the hyperboloids asymptotically approach to. The real axes are the radii of the
    hyperboloids bases, i.e. the radii of the cross sections of the hyperboloids with the xy plane.
    This geometry contains two more parameters, the imaginary radius \f$c\f$ and the radius of the
    hyperboloid top surface, but these values are determined by the other input parameters. */
class HyperboloidShellGeometry : public AxGeometry
{
    ITEM_CONCRETE(HyperboloidShellGeometry, AxGeometry, "a hyperboloid shell geometry")
        ATTRIBUTE_TYPE_DISPLAYED_IF(HyperboloidShellGeometry, "Level2")

        PROPERTY_DOUBLE(outerRadialExtent, "the radial extent of the hyperboloid outer wall")
        ATTRIBUTE_QUANTITY(outerRadialExtent, "length")
        ATTRIBUTE_MIN_VALUE(outerRadialExtent, "]0")

        PROPERTY_DOUBLE(outerOpeningAngle, "the half opening angle of the hyperboloid outer wall")
        ATTRIBUTE_QUANTITY(outerOpeningAngle, "posangle")
        ATTRIBUTE_MIN_VALUE(outerOpeningAngle, "]0 deg")
        ATTRIBUTE_MAX_VALUE(outerOpeningAngle, "90 deg[")

        PROPERTY_DOUBLE(outerRealAxis, "the real axis of the hyperboloid outer wall")
        ATTRIBUTE_QUANTITY(outerRealAxis, "length")
        ATTRIBUTE_MIN_VALUE(outerRealAxis, "]0")

        PROPERTY_DOUBLE(innerOpeningAngle, "the half opening angle of the hyperboloid inner wall")
        ATTRIBUTE_QUANTITY(innerOpeningAngle, "posangle")
        ATTRIBUTE_MIN_VALUE(innerOpeningAngle, "]0 deg")
        ATTRIBUTE_MAX_VALUE(innerOpeningAngle, "90 deg[")

        PROPERTY_DOUBLE(innerRealAxis, "the real axis of the hyperboloid inner wall")
        ATTRIBUTE_QUANTITY(innerRealAxis, "length")
        ATTRIBUTE_MIN_VALUE(innerRealAxis, "]0")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function calculates some frequently used values, including the normalization factor
        \f$A\f$, which is set by the normalization condition that total mass equals one, resulting
        in \f[ A=\frac{3}{2\pi}\frac{1}{z_{\text{max}}
        (2a_{\text{out}}^{2}+b_{\text{out}}^{2}-2a_{\text{in}}^{2}-b_{\text{in}}^{2})}, \f] where
        \f$b\f$ is the radius of the hyperboloid top surface (i.e. cross section with \f$
        z_{\text{max}} \f$ plane), determined by the radial extent and the half opening angle. */
    void setupSelfBefore() override;

    //======================== Other Functions =======================

public:
    /** This function returns the density \f$\rho\f$ at the cylindrical radius \f$R\f$ and height
        \f$z\f$. For the present geometry, it returns a constant density for positions inside the
        hyperboloid shell, and zero for positions outside the hyperboloid shell. */
    double density(double R, double z) const override;

    /** This function generates a random position from the probability density \f$p({\bf{r}})\,
        {\text{d}}{\bf{r}} = \rho({\bf{r}})\, {\text{d}}{\bf{r}}\f$. Since the current
        implementation allows only a constant density distribution, a random position can easily be
        constructed. The function repeatedly generates a uniform random position in the cylinder
        enveloping the geometry, and accepts the first position that happens to fall inside the
        nonzero-density region of the actual geometry. */
    Position generatePosition() const override;

    /** This function returns the radial surface density, i.e. the integration of the density along
        a line in the equatorial plane starting at the centre of the coordinate system, \f[
        \Sigma_R = \int_0^\infty \rho(R,0,0)\,{\text{d}}r. \f] In the case of hyperboloid shell
        geometry with constant density this is simply \f[ \Sigma_R =
        A(a_{\text{max}}-a_{\text{min}}). \f] */
    double SigmaR() const override;

    /** This function returns the Z-axis surface density, i.e. the integration of the density along
        the entire Z-axis, \f[ \Sigma_Z = \int_{-\infty}^\infty \rho(0,0,z)\, {\text{d}}z. \f] For
        the hyperboloid shell geometry this integral is simply zero. */
    double SigmaZ() const override;

    //======================== Data Members ========================

private:
    // aliases to discoverable data members for ease of notation and backwards compatibility
    const double& _Dout{_outerRadialExtent};
    const double& _DeltaOut{_outerOpeningAngle};
    const double& _aout{_outerRealAxis};
    const double& _DeltaIn{_innerOpeningAngle};
    const double& _ain{_innerRealAxis};

    // data members initialized during setup
    double _zmax{0.};
    double _bout{0.};
    double _cout{0.};
    double _Din{0.};
    double _bin{0.};
    double _cin{0.};
    double _A{0.};
};

////////////////////////////////////////////////////////////////////

#endif
