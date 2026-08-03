/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef PARABOLOIDGEOMETRY_HPP
#define PARABOLOIDGEOMETRY_HPP

#include "AxGeometry.hpp"

////////////////////////////////////////////////////////////////////

/** The ParaboloidGeometry class is a subclass of the AxGeometry class and describes the geometry
    defined by a double paraboloid surface. It may be used to represent dusty outflows.

    \image html ParaboloidGeometry.png

    The current implementation allows only a constant density distribution bounded by the
    paraboloid surface \f[ z=\frac{\rho^2}{a^2} \f] and two horizontal planes parallel to the xy
    plane, \f$\pm (z_{\text{max}}+z_{0})\f$. There are three free parameters describing this dust
    geometry: the radial extent \f$D\f$, the half opening angle \f$\Delta\f$, and the offset of the
    paraboloid vertices in the z direction \f$z_{0}\f$. The radial extent corresponds to the length
    of the line connecting the paraboloid vertex and its top for a given half opening angle of the
    paraboloid; the half opening angle is measured from the z axis to the radial extent line. The
    paraboloid height, radius and curvature level are determined by the radial extent and half
    opening angle. */
class ParaboloidGeometry : public AxGeometry
{
    ITEM_CONCRETE(ParaboloidGeometry, AxGeometry, "a paraboloid geometry")
        ATTRIBUTE_TYPE_DISPLAYED_IF(ParaboloidGeometry, "Level2")

        PROPERTY_DOUBLE(radialExtent, "the radial extent of the paraboloid")
        ATTRIBUTE_QUANTITY(radialExtent, "length")
        ATTRIBUTE_MIN_VALUE(radialExtent, "]0")

        PROPERTY_DOUBLE(openingAngle, "the half opening angle of the paraboloid")
        ATTRIBUTE_QUANTITY(openingAngle, "posangle")
        ATTRIBUTE_MIN_VALUE(openingAngle, "]0 deg")
        ATTRIBUTE_MAX_VALUE(openingAngle, "90 deg[")

        PROPERTY_DOUBLE(offsetZ, "the offset of the paraboloid vertices in the z direction")
        ATTRIBUTE_QUANTITY(offsetZ, "length")
        ATTRIBUTE_MIN_VALUE(offsetZ, "[0")
        ATTRIBUTE_DEFAULT_VALUE(offsetZ, "0")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function calculates some frequently used values, including the normalization factor
        \f$A\f$, which is set by the normalization condition that total mass equals one, resulting
        in \f[ A=\frac{1}{\pi*R_{p}^{2}*z_{max}}. \f] */
    void setupSelfBefore() override;

    //======================== Other Functions =======================

public:
    /** This function returns the density \f$\rho\f$ at the cylindrical radius \f$R\f$ and height
        \f$z\f$. For the present geometry, it returns a constant density for positions inside the
        paraboloid, and zero for positions outside the paraboloid. */
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
        z_{max}. \f] */
    double SigmaZ() const override;

    //======================== Data Members ========================

private:
    // aliases to discoverable data members for ease of notation and backwards compatibility
    const double& _D{_radialExtent};
    const double& _Delta{_openingAngle};
    const double& _z0{_offsetZ};

    // data members initialized during setup
    double _Rp{0.};
    double _zmax{0.};
    double _a{0.};
    double _A{0.};
};

////////////////////////////////////////////////////////////////////

#endif
