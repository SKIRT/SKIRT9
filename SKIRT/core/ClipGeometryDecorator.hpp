/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef CLIPGEOMETRYDECORATOR_HPP
#define CLIPGEOMETRYDECORATOR_HPP

#include "Geometry.hpp"

////////////////////////////////////////////////////////////////////

/** The abstract ClipGeometryDecorator class implements a decorator that adjusts another geometry
    by setting the density equal to zero inside or outside a region defined in a subclass. Each
    ClipGeometryDecorator subclass must implement the virtual functions dimension() and inside().
    The decorator increases the density in the remaining region with a constant factor to ensure
    that the total mass remains equal to one. The current implementation does not properly adjust
    the surface densities along the coordinate axes for the mass taken away by the cavity. */
class ClipGeometryDecorator : public Geometry
{
    /** The enumeration type indicating which region to remove (Inside or Outside). */
    ENUM_DEF(Remove, Inside, Outside)
        ENUM_VAL(Remove, Inside, "the inner region (creating a cavity)")
        ENUM_VAL(Remove, Outside, "the outer region (cropping)")
    ENUM_END()

    ITEM_ABSTRACT(ClipGeometryDecorator, Geometry, "a decorator that clips another geometry")
        ATTRIBUTE_TYPE_DISPLAYED_IF(ClipGeometryDecorator, "Level2")

        PROPERTY_ITEM(geometry, Geometry, "the geometry to be clipped")

        PROPERTY_ENUM(remove, Remove, "the region to be removed")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function estimates the fraction \f$\chi\f$ of the mass from the original model taken
        away by the clipping. It samples the density of the geometry being decorated, and counts
        the number of generated positions that fall in the removed region. This value is used to
        renormalize the decorated density distribution to unity: the factor by which the original
        density has to be multiplied is simply \f$1/(1-\chi)\f$. */
    void setupSelfAfter() override;

    //======================== Other Functions =======================

public:
    /** This function returns the density \f$\rho({\bf{r}})\f$ at the position \f${\bf{r}}\f$. It
        is zero in the removed region, and equal to the density of the geometry being decorated
        elsewhere, after an adjustment is made to account for the clipping. */
    double density(Position bfr) const override;

    /** This function generates a random position from the geometry, by drawing a random point from
        the three-dimensional probability density \f$p({\bf{r}})\, {\text{d}}{\bf{r}} =
        \rho({\bf{r}})\, {\text{d}}{\bf{r}}\f$. It repeatedly calls the density() function for the
        geometry being decorated until a position is returned that does not lie in the removed
        region. */
    Position generatePosition() const override;

    /** This function returns the X-axis surface density, i.e. the integration of the density along
        the entire X-axis, \f[ \Sigma_X = \int_{-\infty}^\infty \rho(x,0,0)\,{\text{d}}x. \f] It
        returns the corresponding value of the geometry being decorated after re-normalization. */
    double SigmaX() const override;

    /** This function returns the Y-axis surface density, i.e. the integration of the density along
        the entire Y-axis, \f[ \Sigma_Y = \int_{-\infty}^\infty \rho(0,y,0)\,{\text{d}}y. \f] It
        returns the corresponding value of the geometry being decorated after re-normalization. */
    double SigmaY() const override;

    /** This function returns the Z-axis surface density, i.e. the integration of the density along
        the entire Z-axis, \f[ \Sigma_Z = \int_{-\infty}^\infty \rho(0,0,z)\,{\text{d}}z. \f] It
        returns the corresponding value of the geometry being decorated after re-normalization. */
    double SigmaZ() const override;

protected:
    /** This pure virtual function, to be implemented by a subclass, returns true if the specified
        position is inside the boundary defined by the subclass, i.e. the point is in the region
        that would be carved away when creating a cavity, or in the region that would be retained
        when cropping. */
    virtual bool inside(Position bfr) const = 0;

    /** This function returns the normalization factor calculated by this class during setup. */
    double norm() const;

    //======================== Data Members ========================

private:
    // values initialized during setup
    double _norm{0.};
};

////////////////////////////////////////////////////////////////////

#endif
