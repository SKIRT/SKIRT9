/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef CUBICSPLINESMOOTHINGKERNEL_HPP
#define CUBICSPLINESMOOTHINGKERNEL_HPP

#include "SmoothingKernel.hpp"
#include "Array.hpp"

////////////////////////////////////////////////////////////////////

/** The CubicSpineSmoothingKernel class is a subclass of the abstract SmoothingKernel class, and
    describes smoothing kernels defined by the standard cubic spline density, \f[ W(u) =
    \begin{cases}\; \dfrac{8\,(1-6\,u^2+6\,u^3)}{\pi} & \quad\text{if }0\leq u \leq \tfrac12, \\ \;
    \dfrac{16\,(1-u)^3}{\pi} & \quad \text{if }\tfrac12 \leq u \leq 1, \\ \; 0 & \quad \text{else}.
    \end{cases} \f] It can be checked that this function is continuous at \f$u=\tfrac12\f$ and that
    it satisfies the required normalization \f[ 4\pi \int_0^\infty W(u)\, u^2\, {\text{d}}u = 1.
    \f] */
class CubicSplineSmoothingKernel : public SmoothingKernel
{
    ITEM_CONCRETE(CubicSplineSmoothingKernel, SmoothingKernel, "a cubic spline smoothing kernel")
    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function sets up a grid that will be used to sample random radii from the smoothing
        kernel. */
    void setupSelfBefore() override;

    //======================== Other Functions =======================

public:
    /** This function returns the density \f$W(u)\f$ of the smoothing kernel as a function of the
        normalized radius \f$u\f$. Just implements the analytical formula. */
    double density(double u) const override;

    /** This function generates a random normalized radius \f$u\f$ from the smoothing kernel, by
        drawing a number from the one-dimensional probability density \f$p(u)\,{\text{d}}u =
        4\pi\,W(u)\,u^2\, {\text{d}}u\f$. This is accomplished by generating a uniform deviate
        \f${\cal{X}}\f$, and solving the equation \f[ {\cal{X}} = 4\pi\,W(u)\,u^2\, {\text{d}}u \f]
        for \f$u\f$. For the cubic spline smoothing kernel, we use a precomputed grid with values
        on which we interpolate to solve this equation. */
    double generateRadius() const override;

    //======================== Data Members ========================

private:
    // data members initialized during setup
    Array _Xv;
};

////////////////////////////////////////////////////////////////////

#endif
