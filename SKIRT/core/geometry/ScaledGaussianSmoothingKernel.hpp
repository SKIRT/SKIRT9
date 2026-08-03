/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef SCALEDGAUSSIANSMOOTHINGKERNEL_HPP
#define SCALEDGAUSSIANSMOOTHINGKERNEL_HPP

#include "Array.hpp"
#include "SmoothingKernel.hpp"

////////////////////////////////////////////////////////////////////

/** An instance of the ScaledGaussianSmoothingKernel describes a scaled Gaussian smoothing kernel
    with finite support, as presented by Altay and Theuns 2013 (MNRAS 434,748). For the normalized
    radius \f$u=r/h\f$, the kernel profile is given by:

    \f[ W(u) = \begin{cases}\; \mathcal{N}\,\exp(-\frac{u^2}{2\sigma^2}) & \quad\text{if }0\leq u
    \leq 1, \\ \; 0 & \quad \text{else}. \end{cases} \f]

    where

    \f[ \begin{aligned} \sigma &= \frac{1}{2\sqrt{2}\,\pi^{1/6}} \\ \mathcal{N} &=
    \frac{8}{\pi}\,\left[\mathrm{erf}(t)-\frac{2t\exp(-t^2)}{\sqrt{\pi}} \right]^{-1} \quad
    \text{with} \quad t=2\pi^{1/6}. \end{aligned} \f]

    It can be verified that this function satisfies the required normalization \f[ 4\pi
    \int_0^\infty W(u)\, u^2\, {\text{d}}u = 1. \f]

    With this scaling and cutoff, the Gaussian profile approximates the standard cubic spline
    kernel profile to within about three percent for all radii. The reason for using a Gaussion
    kernel instead of the standard cubic spline kernel in some applications is that a spherical
    Gaussian profile (assuming infinite support, i.e. not cut off at the smoothing length) can be
    separated into Gaussian component profiles along each of the coordinate axes, It thus becomes
    easy to calculate the mass inside a cuboidal box (such as a grid cell) or to determine the
    surface density for the projection on a rectangle (such as a pixel). */
class ScaledGaussianSmoothingKernel : public SmoothingKernel
{
    ITEM_CONCRETE(ScaledGaussianSmoothingKernel, SmoothingKernel, "a scaled Gaussian smoothing kernel")
    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function sets up a grid that will be used to sample random radii from the smoothing
        kernel. */
    void setupSelfBefore() override;

    //======================== Other Functions =======================

public:
    /** This function returns the density \f$W(u)\f$ of the smoothing kernel as a function of the
        normalized radius \f$u\f$. It just implements the analytical formula given in the class
        header. */
    double density(double u) const override;

    /** This function returns the column density \f$\Sigma(q) = 2 \int_{q}^1 \frac{W(u)\,u
        \,{\text{d}}u} {\sqrt{u^2-q^2}}\f$ of the smoothing kernel as a function of the normalized
        impact radius \f$q=r_\text{i}/h\f$. For the scaled Gaussian smoothing kernel, we obtain
        \f[\Sigma(q) = \begin{cases} \; {\cal{N}}\,\sqrt{2\pi}\,\sigma
        \exp\left(-\frac{q^2}{2\sigma^2}\right) \mathrm{erf}\left(\frac{\sqrt{1-q^2}}
        {\sqrt{2}\sigma}\right) & \quad{\text{if }} 0\leq q\leq 1, \\ \; 0 & \quad{\text{else}}.
        \end{cases} \f] */
    double columnDensity(double q) const override;

    /** This function generates a random normalized radius \f$u\f$ from the smoothing kernel, by
        drawing a number from the one-dimensional probability density \f$p(u)\,{\text{d}}u =
        4\pi\,W(u)\,u^2\, {\text{d}}u\f$. This is accomplished by generating a uniform deviate
        \f${\cal{X}}\f$, and solving the equation \f[ {\cal{X}} = \int_0^u 4\pi\,W(u')\,u'^2\,
        {\text{d}}u' \f] for \f$u\f$. For the scaled gaussian smoothing kernel, we use a
        precomputed grid with values on which we interpolate to solve this equation. */
    double generateRadius() const override;

    //======================== Data Members ========================

private:
    // data members initialized during setup
    Array _Xv;
};

////////////////////////////////////////////////////////////////////

#endif
