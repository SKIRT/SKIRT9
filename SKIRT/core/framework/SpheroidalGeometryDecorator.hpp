/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef SPHEROIDALGEOMETRYDECORATOR_HPP
#define SPHEROIDALGEOMETRYDECORATOR_HPP

#include "AxGeometry.hpp"
#include "SpheGeometry.hpp"

////////////////////////////////////////////////////////////////////

/** The SpheroidalGeometryDecorator class is a geometry decorator that constructs a spheroidal
    geometry based on a spherical geometry. The properties of a SpheroidalGeometryDecorator object
    are a reference to the SpheGeometry object being decorated and the flattening parameter
    \f$q\f$. If the original spherical geometry is characterized by the density profile \f$
    \rho_{\text{s}}(r) \f$, the new geometry has as density \f[ \rho(R,z) = \frac{1}{q}\,
    \rho_{\text{s}}\left(\sqrt{R^2 + \frac{z^2}{q^2}}\right). \f] This new geometry is also
    normalized to one.

    Note that the flattening parameter can have any value \f$q>0\f$. For \f$0<q<1\f$, the geometry
    is wider than it is high. For \f$q>1\f$, the geometry is higher than it is wide. */
class SpheroidalGeometryDecorator : public AxGeometry
{
    ITEM_CONCRETE(SpheroidalGeometryDecorator, AxGeometry,
                  "a decorator that constructs a spheroidal variant of any spherical geometry")

        PROPERTY_ITEM(geometry, SpheGeometry, "the spherical geometry to be made spheroidal")

        PROPERTY_DOUBLE(flattening, "the flattening parameter q")
        ATTRIBUTE_MIN_VALUE(flattening, "]0")
        ATTRIBUTE_DEFAULT_VALUE(flattening, "1")

    ITEM_END()

    //======================== Other Functions =======================

public:
    /** This function returns the density \f$\rho(R,z)\f$ at the cylindrical radius \f$R\f$ and
        the height \f$z\f$. It just implements the analytical formula. */
    double density(double R, double z) const override;

    /** This function generates a random position from the geometry, by drawing a random
        point from the three-dimensional probability density \f$p({\bf{r}})\, {\text{d}}{\bf{r}} =
        \rho({\bf{r}})\, {\text{d}}{\bf{r}}\f$. It first generates a random position
        \f${\bf{r}}_{\text{s}}\f$ by calling the generatePosition() function of the geometry
        being decorated and applies a simple linear transformation to the coordinates, \f$x = x_{\text{s}},
        y = y_{\text{s}}, z = q\,z_{\text{s}}\f$. */
    Position generatePosition() const override;

    /** This function returns the radial surface density, i.e. the integration of
        the density along a line in the equatorial plane starting at the centre of the coordinate
        system, \f[ \Sigma_R = \int_0^\infty \rho(R,0)\,{\text{d}}R. \f] We easily obtain
        \f[ \Sigma_R = \frac{1}{q} \int_0^\infty \rho_{\text{orig}}(R)\,{\text{d}}R =
        \frac{1}{q}\,\Sigma_{r,{\text{orig}}}. \f] */
    double SigmaR() const override;

    /** This function returns the Z-axis surface density, i.e. the integration of the density
        along the entire Z-axis, \f[ \Sigma_Z = \int_{-\infty}^\infty \rho(0,0,z)\,{\text{d}}z. \f]
        We easily obtain \f[ \Sigma_Z = \frac{2}{q} \int_0^\infty \rho_{\text{orig}}
        \left(\frac{z}{q}\right)\,{\text{d}}z =
        2 \int_0^\infty \rho_{\text{orig}}(r)\,{\text{d}}r = 2\,\Sigma_{r,{\text{orig}}}. \f] */
    double SigmaZ() const override;

    //======================== Data Members ========================

private:
    // alias to discoverable data member for ease of notation and backwards compatibility
    const double& _q{_flattening};
};

////////////////////////////////////////////////////////////////////

#endif
