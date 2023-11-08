
/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef REDISTRIBUTEGEOMETRYDECORATOR_HPP
#define REDISTRIBUTEGEOMETRYDECORATOR_HPP

#include "Geometry.hpp"

////////////////////////////////////////////////////////////////////

/** The abstract RedistributeGeometryDecorator class implements a decorator that adjusts another
    geometry by multiplying the density with some abstract weight function \f[ \rho'({\bf{r}}) = n
    \rho({\bf{r}}) w({\bf{r}}). \f] Each RedistributeGeometryDecorator subclass must implement the
    virtual functions dimension(), weight() and maxWeight(). The decorator renormalizes the
    distribution by using the importance sampling method from the original distribution. The random
    positions are generated using the rejection method with the reference distribution the original
    density. The current implementation does not properly adjust the surface densities along the
    coordinate axes. */
class RedistributeGeometryDecorator : public Geometry
{
    ITEM_ABSTRACT(RedistributeGeometryDecorator, Geometry, "a decorator that redistributes another geometry")
        ATTRIBUTE_TYPE_DISPLAYED_IF(RedistributeGeometryDecorator, "Level2")

        PROPERTY_ITEM(geometry, Geometry, "the geometry for which the density will be redistributed")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function calculates the norm using the importance sampling method from the original
        distribution. It also calculates the factor \f$\frac{c}{n}\f$ for the rejection method,
        which is just the maxWeight(). */
    void setupSelfAfter() override;

    //======================== Other Functions =======================

public:
    /** This function returns the normalized density \f$n\rho({\bf{r}})w({\bf{r}})\f$ at the
        position \f${\bf{r}}\f$. */
    double density(Position bfr) const override;

    /** This function generates a random position from the geometry using the rejection method. It
        draws a random point from the probability density \f$\rho({\bf{r}})\f$, this point is
        accepted if \f$t=\xi \frac{\max_{\bf{r}}{(w({\bf{r}}))}}{w({\bf{r}})}\le 1.\f$ */
    Position generatePosition() const override;

    /** This function returns the X-axis surface density, i.e. the integration of the density along
        the entire X-axis, \f[ \Sigma_X = \int_{-\infty}^\infty \rho(x,0,0)\,{\text{d}}x. \f] For a
        general geometry this decorator will not have an analytical solution for this integral. We
        use the X-axis surface density of the original distribution. */
    double SigmaX() const override;

    /** This function returns the Y-axis surface density, i.e. the integration of the density along
        the entire Y-axis, \f[ \Sigma_Y = \int_{-\infty}^\infty \rho(y,0,0)\,{\text{d}}y. \f] For a
        general geometry this decorator will not have an analytical solution for this integral. We
        use the Y-axis surface density of the original distribution. */
    double SigmaY() const override;

    /** This function returns the Z-axis surface density, i.e. the integration of the density along
        the entire Z-axis, \f[ \Sigma_Z = \int_{-\infty}^\infty \rho(Z,0,0)\,{\text{d}}z. \f] For a
        general geometry this decorator will not have an analytical solution for this integral. We
        use the Z-axis surface density of the original distribution. */
    double SigmaZ() const override;

protected:
    /** This pure virtual function, to be implemented by a subclass, gives the weight function's
        value at a given position. This function does not need to be normalized. */
    virtual double weight(Position bfr) const = 0;

    /** This pure virtual function, to be implemented by a subclass, gives the maximum value of the
        weight function. This value is used in the rejection method. */
    virtual double maxWeight() const = 0;

    //======================== Data Members ========================

private:
    // values initialized during setup
    double _norm{0.};
    double _maxWeight;
};

////////////////////////////////////////////////////////////////////

#endif
