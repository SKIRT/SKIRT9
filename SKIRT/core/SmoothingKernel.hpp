/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef SMOOTHINGKERNEL_HPP
#define SMOOTHINGKERNEL_HPP

#include "SimulationItem.hpp"
#include "StoredTable.hpp"
class Random;
class Box;

//////////////////////////////////////////////////////////////////////

/** SmoothingKernel is an abstract class that represents smoothing kernels that can be used to
    smear out points to smoothed particles. Each smoothing kernel is completely defined by the
    kernel density \f$W(u)\f$ which is just a function of the normalized radius \f$u=r/h\f$ with
    \f$h\f$ the smoothing length. All smoothing kernels are spherically symmetric, have finite
    support, and are normalized to one, i.e. \f[ W(u)=0 \quad \text{for}\; u>1 \qquad \text{and}
    \qquad 4\pi \int_0^1 W(u)\, u^2\, {\text{d}}u = 1. \f] */
class SmoothingKernel : public SimulationItem
{
    ITEM_ABSTRACT(SmoothingKernel, SimulationItem, "a smoothing kernel")
    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function caches the simulation's random generator for use by subclasses and opens
     *  the stored table with the tabulated cumulative kernel. */
    void setupSelfBefore() override;

    //======================== Other Functions =======================

public:
    /** This function returns the density \f$W(u)\f$ of the smoothing kernel as a function of the
        normalized radius \f$u\f$. Subclasses must implement this function appropriately. */
    virtual double density(double u) const = 0;

    /** This function returns the projection of the smoothing kernel along a line with given impact
        radius (smallest distance) to the kernel's center. In formula form, this is the column
        density \f$\Sigma(q) = 2 \int_{q}^1 \frac{W(u)\,u \,{\text{d}}u} {\sqrt{u^2-q^2}}\f$ of the
        smoothing kernel as a function of the normalized impact radius \f$q=r_\text{i}/h\f$, or
        zero for \f$q\ge 1\f$. Subclasses must implement this function appropriately. */
    virtual double columnDensity(double q) const = 0;

    /** This pure virtual function generates a random normalized radius \f$u\f$ from the smoothing
        kernel, by drawing a number from the one-dimensional probability density \f$
        p(u)\,{\text{d}}u = 4\pi\,W(u)\,u^2\, {\text{d}}u \f$. Subclasses must implement this
        function appropriately. */
    virtual double generateRadius() const = 0;

    /** This function returns the portion of the smoothing kernel mass contained in the specified
        axis-aligned bounding box, i.e.: \f[ m_\text{box} = \int_{x_\text{min}}^{x_\text{max}}
        \int_{y_\text{min}}^{y_\text{max}} \int_{z_\text{min}}^{z_\text{max}} \text{d}x \text{d}y
        \text{d}z. \f]

        The implementation in this base class uses the technique described by Baes et al. 2026 in
        prep, reducing the calculation to eight interpolations in the tabulated cumulative kernel
        provided as a stored table resource with the same name as the final derived class. Derived
        classes can override the function to provide another implementation if desired. */
    virtual double massInBox(const Box& box) const;

    /** This function returns true if this smoothing kernel offers the massInBox() feature. The
        implementation in this base class returns true if the cumulative kernel resource file is
        available, and false if not. Derived classes can override the function if desired. */
    virtual bool hasMassInBox() const;

protected:
    /** This function returns the simulation's random generator as a service to subclasses. */
    Random* random() const { return _random; }

    //======================== Data Members ========================

private:
    // data member initialized during setup
    Random* _random{nullptr};

    // the tabulated cumulative kernel used in massInBox()
    StoredTable<3> _cumkernel;
};

//////////////////////////////////////////////////////////////////////

#endif
