/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef UNIFORMSMOOTHINGKERNEL_HPP
#define UNIFORMSMOOTHINGKERNEL_HPP

#include "SmoothingKernel.hpp"

////////////////////////////////////////////////////////////////////

/** The UniformSmoothingKernel class is a subclass of the abstract SmoothingKernel class, and
    describes smoothing kernels with a uniform density and a compact support. They are completely
    described by the simple kernel density, \f[ W(u) = \begin{cases}\;\dfrac{3}{4\pi} &
    \quad{\text{if }} 0 \leq u \leq 1 \\ \; 0 & \quad{\text{else}} \end{cases}. \f] This function
    satisfies the required normalization \f[ 4\pi \int_0^\infty W(u)\, u^2\, {\text{d}}u = 1. \f]
    */
class UniformSmoothingKernel : public SmoothingKernel
{
    ITEM_CONCRETE(UniformSmoothingKernel, SmoothingKernel, "a uniform smoothing kernel")
    ITEM_END()

    //======================== Other Functions =======================

public:
    /** This function returns the density \f$W(u)\f$ of the smoothing kernel as a function of the
        normalized radius \f$u\f$. Just implements the analytical formula. */
    double density(double u) const override;

    /** This function generates a random normalized radius \f$u\f$ from the smoothing kernel, by
        drawing a number from the one-dimensional probability density \f$p(u)\,{\text{d}}u =
        4\pi\,W(u)\,u^2\, {\text{d}}u\f$. This is accomplished by generating a uniform deviate
        \f${\cal{X}}\f$, and solving the equation \f[ {\cal{X}} = 4\pi\,W(u)\,u^2\, {\text{d}}u \f]
        for \f$u\f$. For the uniform smoothing kernel, we obtain the simple expression \f$ u =
        \sqrt[3]{\cal{X}} \f$. */
    double generateRadius() const override;
};

////////////////////////////////////////////////////////////////////

#endif
