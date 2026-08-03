/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef UNIFORMSMOOTHINGKERNEL_HPP
#define UNIFORMSMOOTHINGKERNEL_HPP

#include "SmoothingKernel.hpp"

////////////////////////////////////////////////////////////////////

/** The UniformSmoothingKernel class is a subclass of the abstract SmoothingKernel class, and
    describes a smoothing kernel with a uniform density and a compact support. It is completely
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
        normalized radius \f$u\f$. It just implements the analytical formula given in the class
        header. */
    double density(double u) const override;

    /** This function returns the column density \f$\Sigma(q) = 2 \int_{q}^1 \frac{W(u)\,u
        \,{\text{d}}u} {\sqrt{u^2-q^2}}\f$ of the smoothing kernel as a function of the normalized
        impact radius \f$q=r_\text{i}/h\f$. For the uniform smoothing kernel, we obtain
        \f[\Sigma(q) = \begin{cases} \;\dfrac{3}{2\pi} \sqrt{1-q^2} & \quad\text{if } 0 \leq q \leq
        1, \\ \; 0 & \quad{\text{else}}. \end{cases} \f] */
    double columnDensity(double q) const override;

    /** This function generates a random normalized radius \f$u\f$ from the smoothing kernel by
        drawing a number from the one-dimensional probability density \f$p(u)\,{\text{d}}u =
        4\pi\,W(u)\,u^2\, {\text{d}}u\f$. This is accomplished by generating a uniform deviate
        \f${\cal{X}}\f$, and solving the equation \f[ {\cal{X}} = 4\pi\,W(u)\,u^2\, {\text{d}}u \f]
        for \f$u\f$. For the uniform smoothing kernel, we obtain the simple expression \f$ u =
        \sqrt[3]{\cal{X}} \f$. */
    double generateRadius() const override;
};

////////////////////////////////////////////////////////////////////

#endif
