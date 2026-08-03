/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef COMBINEGEOMETRYDECORATOR_HPP
#define COMBINEGEOMETRYDECORATOR_HPP

#include "Geometry.hpp"

//////////////////////////////////////////////////////////////////////

/** The CombineGeometryDecorator class implements a decorator that combines two arbitrary
    geometries (with a relative weight factor) into an object that behaves as a single geometry. */
class CombineGeometryDecorator : public Geometry
{
    ITEM_CONCRETE(CombineGeometryDecorator, Geometry, "a decorator that combines two different geometries")

        PROPERTY_ITEM(firstGeometry, Geometry, "the first geometry")

        PROPERTY_DOUBLE(firstWeight, "the weight of the first geometry")
        ATTRIBUTE_MIN_VALUE(firstWeight, "]0")
        ATTRIBUTE_DEFAULT_VALUE(firstWeight, "1")

        PROPERTY_ITEM(secondGeometry, Geometry, "the second geometry")

        PROPERTY_DOUBLE(secondWeight, "the weight of the second geometry")
        ATTRIBUTE_MIN_VALUE(secondWeight, "]0")
        ATTRIBUTE_DEFAULT_VALUE(secondWeight, "1")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function normalizes the relative weights of the two geometries. */
    void setupSelfBefore() override;

    //======================== Other Functions =======================

public:
    /** This function returns the dimension of the combined geometry, which depends on the (lack of)
        symmetry in the geometries of its stellar components. A value of 1 means spherical
        symmetry, 2 means axial symmetry and 3 means none of these symmetries. The geometry
        with the least symmetry (i.e. the highest dimension) determines the result for
        the whole system. */
    int dimension() const override;

    /** This function returns the density \f$\rho({\bf{r}})\f$ at the position
        \f${\bf{r}}\f$. It adds the contribution of the two components, weighted by their
        weights. */
    double density(Position bfr) const override;

    /** This function generates a random position from the geometry, by drawing a random
        point from the three-dimensional probability density \f$p({\bf{r}})\, {\text{d}}{\bf{r}} =
        \rho({\bf{r}})\, {\text{d}}{\bf{r}}\f$. It first generates a random component, and
        subsequently generates a random position from this component by calling the corresponding
        generatePosition() function. */
    Position generatePosition() const override;

    /** This function returns the X-axis surface density, i.e. the integration of the density
        along the entire X-axis, \f[ \Sigma_X = \int_{-\infty}^\infty \rho(x,0,0)\,{\text{d}}x. \f]
        It adds the contribution of the two components, weighted by their weights. */
    double SigmaX() const override;

    /** This function returns the Y-axis surface density, i.e. the integration of the density
        along the entire Y-axis, \f[ \Sigma_Y = \int_{-\infty}^\infty \rho(0,y,0)\,{\text{d}}y. \f]
        It adds the contribution of the two components, weighted by their weights. */
    double SigmaY() const override;

    /** This function returns the Z-axis surface density, i.e. the integration of the density
        along the entire Z-axis, \f[ \Sigma_Z = \int_{-\infty}^\infty \rho(0,0,z)\,{\text{d}}z. \f]
        It adds the contribution of the two components, weighted by their weights. */
    double SigmaZ() const override;

    //======================== Data Members ========================

private:
    double _w1{0}, _w2{0};  // the normalized weights
};

////////////////////////////////////////////////////////////////

#endif
