/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef GEOMETRY_HPP
#define GEOMETRY_HPP

#include "Position.hpp"
#include "SimulationItem.hpp"
class Random;

//////////////////////////////////////////////////////////////////////

/** Geometry is the general abstract class that describes the density distribution of a primary
    source or of some obscuring medium. It is an abstract class, which essentially describes the
    interface that needs to be implemented by any subclass. There are two key member functions that
    each Geometry subclass should have: a function returning the density \f$\rho({\bf{r}})\f$ and a
    function drawing a random position from this geometry. We normalize each geometry such that the
    total mass is equal to one, i.e. \f[ \iiint \rho({\bf{r}})\, {\text{d}}{\bf{r}} = 1.\f] */
class Geometry : public SimulationItem
{
    ITEM_ABSTRACT(Geometry, SimulationItem, "a geometry")
    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function caches the simulation's random generator for use by subclasses. */
    void setupSelfBefore() override;

    //======================== Other Functions =======================

public:
    /** This pure virtual function returns the dimension of the geometry, which depends on its
        (lack of) symmetry. A value of 1 means spherical symmetry, 2 means axial symmetry and 3
        means none of these symmetries. The function's implementation must be provided in a
        subclass. */
    virtual int dimension() const = 0;

    /** This pure virtual function returns the density \f$\rho({\bf{r}})\f$ at the position
        \f${\bf{r}}\f$. */
    virtual double density(Position bfr) const = 0;

    /** This pure virtual function generates a random position from the geometry, by
        drawing a random point from the three-dimensional probability density \f$p({\bf{r}})\,
        {\text{d}}{\bf{r}} = \rho({\bf{r}})\, {\text{d}}{\bf{r}}\f$. */
    virtual Position generatePosition() const = 0;

    /** This pure virtual function returns the X-axis surface density, i.e. the integration of
        the density along the entire X-axis, \f[ \Sigma_X = \int_{-\infty}^\infty \rho(x,0,0)\,
        {\text{d}}x. \f] */
    virtual double SigmaX() const = 0;

    /** This pure virtual function returns the Y-axis surface density, i.e. the integration of
        the density along the entire Y-axis, \f[ \Sigma_Y = \int_{-\infty}^\infty \rho(0,y,0)\,
        {\text{d}}y. \f] */
    virtual double SigmaY() const = 0;

    /** This pure virtual function returns the Z-axis surface density, i.e. the integration of
        the density along the entire Z-axis, \f[ \Sigma_Z = \int_{-\infty}^\infty \rho(0,0,z)\,
        {\text{d}}z. \f] */
    virtual double SigmaZ() const = 0;

    //======================== Other Functions =======================

protected:
    /** This function returns the simulation's random generator as a service to subclasses. */
    Random* random() const { return _random; }

    //======================== Data Members ========================

private:
    // data member initialized during setup
    Random* _random{nullptr};
};

////////////////////////////////////////////////////////////////////

#endif
