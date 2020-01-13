/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef HYPERBOLOIDGEOMETRY_HPP
#define HYPERBOLOIDGEOMETRY_HPP

#include "AxGeometry.hpp"

////////////////////////////////////////////////////////////////////

/** The HyperboloidGeometry class is a subclass of the AxGeometry class and describes the geometry
    defined by a hyperboloid surface. It may be used to represent dusty outflows, extending
    vertically and then flaring towards the outside, asymptotically approaching a cone surface.

    \image html HyperboloidGeometry.png

    The current implementation allows only a constant density distribution bounded by the
    hyperboloid surface \f[ z=\frac{c}{a}\sqrt{\rho^2-a^2} \f] and two horizontal planes parallel
    to the xy plane, \f$\pm z_{\text{max}}\f$. There are four free parameters describing this dust
    geometry: the radial extent \f$D\f$, the half opening angle \f$\Delta\f$, the real axis \f$a\f$
    and the minimum radius \f$r_{\text{min}}\f$. The first two correspond to the length of the side
    and half opening angle of the cone to which the hyperboloid asymptotically approaches. The real
    axis is the radius of the hyperboloid base, i.e. the radius of the cross section of the
    hyperboloid with the xy plane. The minimum radius defines a spherical cavity centered on the
    coordinate origin. The hyperboloid formula contains two more parameters, the imaginary radius
    \f$c\f$ and the radius of the hyperboloid top surface, but these values are determined by the
    other input parameters.

    If the dusty system under consideration is surrounding an AGN or another source which is
    luminous enough to heat the dust to sublimation temperature, the inner radius should correspond
    to the sublimation radius which scale as \f$ r_{\text{min}} \propto L(\theta)^{0.5}\f$
    (Barvainis, 1987, ApJ, 320, 537, eq (5)). If the primary source assumes anisotropic emission,
    the inner radius should follow the same dependence as the distribution of the primary source
    luminosity. Otherwise, dust temperature on the inner boundary of the geometry is very likely to
    be under- or over-estimated. Thus, if the NetzerAccretionDiskGeometry distribution is chosen to
    describe the primary source emission, it is recommended to enable the "reshape inner radius"
    option. The inner radius will then be set by the following formula: \f[ r_{\text{min}} \propto
    (\cos\theta\,(2\cos\theta+1))^{0.5}.\f] This should allow dust to approach all the way to the
    primary central source in the equatorial plane. However, due to the finite resolution of dust
    cells, it may happen that some of the innermost cells end up with unphysically high
    temperatures. For this reason, there is an additional input parameter, the cutoff radius
    \f$r_{\text{cut}}\f$. The value of the cutoff radius is usually found after a few
    trial-and-error experiments by inspecting temperature distribution maps, until the inner wall
    of the geometry is at the expected sublimation temperature for a given dust population.

    The total dust mass of the model corresponds to the mass of the original geometry, before the
    inner wall is reshaped to account for anisotropy. The deviation may or may not be significant,
    depending on the volume of the entire hyperboloid compared to the central cavity; it is the
    user's responsibility to check if it is within satisfactory limits. */
class HyperboloidGeometry : public AxGeometry
{
    ITEM_CONCRETE(HyperboloidGeometry, AxGeometry, "a hyperboloid geometry")
        ATTRIBUTE_TYPE_DISPLAYED_IF(HyperboloidGeometry, "Level2")

        PROPERTY_DOUBLE(radialExtent, "the radial extent of the hyperboloid")
        ATTRIBUTE_QUANTITY(radialExtent, "length")
        ATTRIBUTE_MIN_VALUE(radialExtent, "]0")

        PROPERTY_DOUBLE(openingAngle, "the half opening angle of the hyperboloid")
        ATTRIBUTE_QUANTITY(openingAngle, "posangle")
        ATTRIBUTE_MIN_VALUE(openingAngle, "]0 deg")
        ATTRIBUTE_MAX_VALUE(openingAngle, "90 deg[")

        PROPERTY_DOUBLE(realAxis, "the real axis of the hyperboloid")
        ATTRIBUTE_QUANTITY(realAxis, "length")
        ATTRIBUTE_MIN_VALUE(realAxis, "]0")

        PROPERTY_DOUBLE(minRadius, "the radius of the central cavity")
        ATTRIBUTE_QUANTITY(minRadius, "length")
        ATTRIBUTE_MIN_VALUE(minRadius, "]0")

        PROPERTY_BOOL(reshapeInnerRadius, "reshape the inner radius according to the Netzer luminosity profile")
        ATTRIBUTE_DEFAULT_VALUE(reshapeInnerRadius, "false")

        PROPERTY_DOUBLE(cutoffRadius, "the inner cutoff radius of the hyperboloid")
        ATTRIBUTE_QUANTITY(cutoffRadius, "length")
        ATTRIBUTE_MIN_VALUE(cutoffRadius, "[0")
        ATTRIBUTE_DEFAULT_VALUE(cutoffRadius, "0")
        ATTRIBUTE_RELEVANT_IF(cutoffRadius, "reshapeInnerRadius")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function calculates some frequently used values, including the normalization factor
        \f$A\f$, which is set by the normalization condition that total mass equals one, resulting
        in \f[ A=\frac{3}{2\pi}\frac{1}{z_{max}(2a^{2}+b^{2}-4 r^{3}_{min})}, \f] where \f$b\f$ is
        the radius of the hyperboloid top surface (i.e. cross section with \f$ z_{\text{max}} \f$
        plane), determined by the radial extent and the half opening angles. */
    void setupSelfBefore() override;

    //======================== Other Functions =======================

public:
    /** This function returns the density \f$\rho\f$ at the cylindrical radius \f$R\f$ and height
        \f$z\f$. For the present geometry, it returns a constant density for positions inside the
        hyperboloid, and zero for positions outside the hyperboloid. */
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
        \Sigma_R = \int_0^\infty \rho(R,0,0)\,{\text{d}}r. \f] In the case of hyperboloid geometry
        with constant density this is simply \f[ \Sigma_R = 2A(a-r_{\text{min}}). \f] */
    double SigmaR() const override;

    /** This function returns the Z-axis surface density, i.e. the integration of the density along
        the entire Z-axis, \f[ \Sigma_Z = \int_{-\infty}^\infty \rho(0,0,z)\, {\text{d}}z. \f] For
        the hyperboloid geometry this integral is simply zero. */
    double SigmaZ() const override;

    //======================== Data Members ========================

private:
    // aliases to discoverable data members for ease of notation and backwards compatibility
    const double& _D{_radialExtent};
    const double& _Delta{_openingAngle};
    const double& _a{_realAxis};
    const double& _rmin{_minRadius};
    const bool& _rani{_reshapeInnerRadius};
    const double& _rcut{_cutoffRadius};

    // data members initialized during setup
    double _zmax{0.};
    double _b{0.};
    double _c{0.};
    double _A{0.};
};

////////////////////////////////////////////////////////////////////

#endif
