/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef PARABOLOIDSHELLGEOMETRY_HPP
#define PARABOLOIDSHELLGEOMETRY_HPP

#include "AxGeometry.hpp"

////////////////////////////////////////////////////////////////////

/** The ParaboloidShellGeometry class is a subclass of the AxGeometry class and describes the
    geometry defined by two paraboloid surfaces. It may be used to represent dusty outflows.

    \image html ParaboloidShellGeometry.png

    The current implementation allows only a constant density distribution, bounded by two
    paraboloid surfaces of the form \f[ z=\frac{\rho^2}{a^2} \f] and two horizontal planes paralel
    to the xy plane, \f$\pm (z_{\text{max}}+z_{\text{in}})\f$. There are five free parameters
    describing this dust geometry: the radial extent of the inner paraboloid \f$D_{\text{in}}\f$,
    the two half opening angles \f$\Delta_{\text{in}}\f$ and \f$\Delta_{\text{out}}\f$, and the two
    offsets of the paraboloid vertices in the z direction \f$z_{\text{in}}\f$ and
    \f$z_{\text{out}}\f$. The radial extent correspond to the length of the line connecting the
    inner paraboloid vertex and its top for a given half opening angle of the paraboloid; the half
    opening angles are measured from the z axis to the radial extent line. The paraboloid height,
    radius and curvature level are determined by the radial extent and half opening angle. */
class ParaboloidShellGeometry : public AxGeometry
{
    ITEM_CONCRETE(ParaboloidShellGeometry, AxGeometry, "a paraboloid shell geometry")
        ATTRIBUTE_TYPE_DISPLAYED_IF(ParaboloidShellGeometry, "Level2")

        PROPERTY_DOUBLE(innerRadialExtent, "the radial extent of the inner paraboloid wall")
        ATTRIBUTE_QUANTITY(innerRadialExtent, "length")
        ATTRIBUTE_MIN_VALUE(innerRadialExtent, "]0")

        PROPERTY_DOUBLE(innerOpeningAngle, "the half opening angle of the inner paraboloid wall")
        ATTRIBUTE_QUANTITY(innerOpeningAngle, "posangle")
        ATTRIBUTE_MIN_VALUE(innerOpeningAngle, "]0 deg")
        ATTRIBUTE_MAX_VALUE(innerOpeningAngle, "90 deg[")

        PROPERTY_DOUBLE(outerOpeningAngle, "the half opening angle of the outer paraboloid wall")
        ATTRIBUTE_QUANTITY(outerOpeningAngle, "posangle")
        ATTRIBUTE_MIN_VALUE(outerOpeningAngle, "]0 deg")
        ATTRIBUTE_MAX_VALUE(outerOpeningAngle, "90 deg[")

        PROPERTY_DOUBLE(innerOffsetZ, "the offset of the inner paraboloid vertices in the z direction")
        ATTRIBUTE_QUANTITY(innerOffsetZ, "length")
        ATTRIBUTE_MIN_VALUE(innerOffsetZ, "[0")
        ATTRIBUTE_DEFAULT_VALUE(innerOffsetZ, "0")

        PROPERTY_DOUBLE(outerOffsetZ, "the offset of the outer paraboloid vertices in the z direction")
        ATTRIBUTE_QUANTITY(outerOffsetZ, "length")
        ATTRIBUTE_MIN_VALUE(outerOffsetZ, "[0")
        ATTRIBUTE_DEFAULT_VALUE(outerOffsetZ, "0")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function calculates some frequently used values, including the normalization factor
        \f$A\f$, which is set by the normalization condition that total mass equals one, resulting
        in \f[ A=\frac{1}{\pi(R_{out}^{2}z_{maxOut}-R_{in}^{2}z_{maxIn})} \f]. */
    void setupSelfBefore() override;

    //======================== Other Functions =======================

public:
    /** This function returns the density \f$\rho\f$ at the cylindrical radius \f$R\f$ and height
        \f$z\f$. For the present geometry, it returns a constant density for positions inside the
        paraboloid shell, and zero for positions outside the paraboloid shell. */
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
        \Sigma_R = \int_0^\infty \rho(R,0,0)\,{\text{d}}r. \f] For the paraboloid geometry this
        integral is simply zero. */
    double SigmaR() const override;

    /** This function returns the Z-axis surface density, i.e. the integration of the density along
        the entire Z-axis, \f[ \Sigma_Z = \int_{-\infty}^\infty \rho(0,0,z)\, {\text{d}}z. \f] In
        the case of paraboloid geometry with constant density this is simply \f[ \Sigma_Z = 2 A
        (z_{\text{in}}-z_{\text{out}}), \f] where \f$z_{\text{in}}\f$ and \f$z_{\text{out}}\f$ are,
        respectively, offsets of the inner and outer paraboloid vertices in the z direction. */
    double SigmaZ() const override;

    //======================== Data Members ========================

private:
    // aliases to discoverable data members for ease of notation and backwards compatibility
    const double& _Din{_innerRadialExtent};
    const double& _DeltaIn{_innerOpeningAngle};
    const double& _DeltaOut{_outerOpeningAngle};
    const double& _zIn{_innerOffsetZ};
    const double& _zOut{_outerOffsetZ};

    // data members initialized during setup
    double _zmaxIn{0.};
    double _Rin{0.};
    double _aIn{0.};
    double _zmaxOut{0.};
    double _Rout{0.};
    double _aOut{0.};
    double _A{0.};
};

////////////////////////////////////////////////////////////////////

#endif
