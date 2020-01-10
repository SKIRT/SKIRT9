/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef TRIAXIALGEOMETRYDECORATOR_HPP
#define TRIAXIALGEOMETRYDECORATOR_HPP

#include "GenGeometry.hpp"
#include "SpheGeometry.hpp"

////////////////////////////////////////////////////////////////////

/** The TriaxialGeometryDecorator class is a geometry decorator that constructs a triaxial geometry
    based on a spherical geometry. The properties of an TriaxialGeometryDecorator object are a
    reference to the SpheGeometry object being decorated and the flattening parameters \f$p\f$ and
    \f$q\f$. If the original spherical geometry is characterized by the density profile \f$
    \rho_{\text{s}}(r) \f$, the new geometry has as density \f[ \rho(x,y,z) = \frac{1}{p\,q}\,
    \rho_{\text{s}}\left(\sqrt{x^2 + \frac{y^2}{p^2} + \frac{z^2}{q^2}}\right). \f] This new
    geometry is also normalized to one. Note that the flattening parameters can have any value
    \f$p>0, q>0\f$. */
class TriaxialGeometryDecorator : public GenGeometry
{
    ITEM_CONCRETE(TriaxialGeometryDecorator, GenGeometry,
                  "a decorator that constructs a triaxial variant of any spherical geometry")

        PROPERTY_ITEM(geometry, SpheGeometry, "the spherical geometry to be made triaxial")

        PROPERTY_DOUBLE(flatteningY, "the flattening parameter p (along the y-axis)")
        ATTRIBUTE_MIN_VALUE(flatteningY, "]0")
        ATTRIBUTE_DEFAULT_VALUE(flatteningY, "1")

        PROPERTY_DOUBLE(flatteningZ, "the flattening parameter q (along the z-axis)")
        ATTRIBUTE_MIN_VALUE(flatteningZ, "]0")
        ATTRIBUTE_DEFAULT_VALUE(flatteningZ, "1")

    ITEM_END()

    //======================== Other Functions =======================

public:
    /** This function returns the density \f$\rho({\bf{r}})\f$ at the position \f${\bf{r}}\f$. It
        just implements the analytical formula. */
    double density(Position bfr) const override;

    /** This function generates a random position from the geometry, by drawing a random
        point from the three-dimensional probability density \f$p({\bf{r}})\, {\text{d}}{\bf{r}} =
        \rho({\bf{r}})\, {\text{d}}{\bf{r}}\f$. It first generates a random position
        \f${\bf{r}}_{\text{s}}\f$ by calling the generatePosition() function of the geometry
        being decorated and applies a simple linear transformation to the coordinates, \f$x = x_{\text{s}},
        y = p\,y_{\text{s}}, z = q\,z_{\text{s}}\f$. */
    Position generatePosition() const override;

    /** This function returns the X-axis surface density, i.e. the integration of the density along
        the entire X-axis, \f[ \Sigma_X = \int_{-\infty}^\infty \rho(x,0,0)\,{\text{d}}x. \f] We
        easily obtain \f[ \Sigma_X = \frac{2}{p\,q} \int_{-\infty}^\infty
        \rho_{\text{orig}}(x)\,{\text{d}}x = \frac{2}{p\,q}\,\Sigma_{r,{\text{orig}}}. \f] */
    double SigmaX() const override;

    /** This function returns the Y-axis surface density, i.e. the integration of the density along
        the entire Y-axis, \f[ \Sigma_Y = \int_{-\infty}^\infty \rho(0,y,0)\,{\text{d}}y. \f] We
        easily obtain \f[ \Sigma_Y = \frac{2}{p\,q} \int_{-\infty}^\infty \rho_{\text{orig}}
        \left(\frac{y}{p}\right)\,{\text{d}}y = \frac{2}{q}\,\Sigma_{r,{\text{orig}}}. \f] */
    double SigmaY() const override;

    /** This function returns the Z-axis surface density, i.e. the integration of the density along
        the entire Z-axis, \f[ \Sigma_Z = \int_{-\infty}^\infty \rho(0,0,z)\,{\text{d}}z. \f] We
        easily obtain \f[ \Sigma_Z = \frac{2}{p\,q} \int_{-\infty}^\infty \rho_{\text{orig}}
        \left(\frac{z}{q}\right)\,{\text{d}}z = \frac{2}{p}\,\Sigma_{r,{\text{orig}}}. \f] */
    double SigmaZ() const override;

    //======================== Data Members ========================

private:
    // aliases to discoverable data members for ease of notation and backwards compatibility
    const double& _p{_flatteningY};
    const double& _q{_flatteningZ};
};

////////////////////////////////////////////////////////////////////

#endif
