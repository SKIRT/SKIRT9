/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef QUARTICSPLINESMOOTHINGKERNEL_HPP
#define QUARTICSPLINESMOOTHINGKERNEL_HPP

#include "Array.hpp"
#include "SmoothingKernel.hpp"

////////////////////////////////////////////////////////////////////

/** The QuarticSplineSmoothingKernel class is a subclass of the abstract SmoothingKernel class, and
    describes a smoothing kernel defined by the quartic spline density, \f[ W(u) =
    \dfrac{15625}{512\pi} \begin{cases} \; 6 u^{4}-\dfrac{12}{5} u^{2}+\dfrac{46}{125} &\quad 0\le
    u \le\dfrac{1}{5} \\ \; -4 u^{4}+8 u^{3}-\dfrac{24}{5} u^{2}+\dfrac{8}{25} u +\dfrac{44}{125}
    &\quad \dfrac{1}{5}\le u \le\dfrac{3}{5} \\ \; u^{4}-4 u^{3}+6 u^{2}-4 u +1 &\quad
    \dfrac{3}{5}\le u \le 1 \\ \; 0 &\quad \text{else}. \end{cases} \f]

    It can be verified that this function is continuous at the break points and that it satisfies
    the required normalization \f[ 4\pi \int_0^\infty W(u)\, u^2\, {\text{d}}u = 1. \f] */
class QuarticSplineSmoothingKernel : public SmoothingKernel
{
    ITEM_CONCRETE(QuarticSplineSmoothingKernel, SmoothingKernel, "a quartic spline smoothing kernel")
    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function sets up a grid that will be used to sample random radii from the smoothing
        kernel. */
    void setupSelfBefore() override;

    //======================== Other Functions =======================

public:
    /** This function returns the density \f$W(u)\f$ of the smoothing kernel as a function of the
        normalized radius \f$u\f$. It implements the formula given in the class header, but using
        floating point coefficients that were precomputed at high precision. */
    double density(double u) const override;

    /** This function returns the column density \f$\Sigma(q) = 2 \int_{q}^1 \frac{W(u)\,u
        \,{\text{d}}u} {\sqrt{u^2-q^2}}\f$ of the smoothing kernel as a function of the normalized
        impact radius \f$q=r_\text{i}/h\f$. For the quartic spline smoothing kernel, this integral
        can be calculated analytically, but the resulting formula is very complicated and becomes
        numerically unstable near some interval borders. Rather than attempting to implement this
        formula, we implement power series approximations over 8 distinct intervals covering the
        complete domain. These series and the corresponding floating point coefficients were
        obtained through a symbolic software package. The approximation has a maximum absolute
        error of \f$10^{-7}\f$, compared to a maximum column density value of \f$\Sigma(q=0)
        \approx 2.3873\f$. */
    double columnDensity(double q) const override;

    /** This function generates a random normalized radius \f$u\f$ from the smoothing kernel. This
        is accomplished by generating a uniform deviate \f${\cal{X}}\f$, and solving the equation
        \f[ {\cal{X}} = \int_0^u 4\pi\,W(u')\,u'^2\, {\text{d}}u' \f] for \f$u\f$. For the quartic
        spline smoothing kernel, this cumulative distribution function can be readily obtained
        through a symbolic software package. We use this result, with precomputed high-precision
        floating point coefficients, to precompute a tabulation of the cumulative distribution
        function on which we interpolate to solve the above equation. */
    double generateRadius() const override;

    //======================== Data Members ========================

private:
    // data members initialized during setup
    Array _Xv;
};

////////////////////////////////////////////////////////////////////

#endif
