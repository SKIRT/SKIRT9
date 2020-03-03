/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef VOIGTPROFILE_HPP
#define VOIGTPROFILE_HPP

#include "Basics.hpp"

////////////////////////////////////////////////////////////////////

/** This namespace offers functions related to the Voigt profile, which is a probability
    distribution given by a convolution of a Cauchy-Lorentz distribution and a Gaussian
    distribution. */
namespace VoigtProfile
{
    /** This function returns (an approximation to) the value of the Voigt function \f$H(a,x)\f$
        defined by \f[ H(a,x) = \frac{a}{\pi} \int_{-\infty}^\infty
        \frac{\mathrm{e}^{-y^2}\,\mathrm{d}y}{(y-x)^2+a^2} \approx \begin{cases} \mathrm{e}^{-x^2}
        & \text{core} \\ \dfrac{a}{\sqrt{\pi}x^2} & \text{wings} \end{cases} \f] where \f$a\f$ is a
        parameter and the function is normalized so that \f$H(a,0) = 1\f$.

        We use the approximation provided by Tasitsiomi 2006 (ApJ, 645, 792-813) and Laursen et al.
        2009 (ApJ, 696, 853-869), which can be written as \f[ H(a,x) = q + \exp(-x^2) \f] where \f[
        q = \begin{cases} 0 & \mathrm{for} \quad z \le 0 \\ \displaystyle
        \frac{a}{\sqrt{\pi}(x^2+1)} \left(1+\frac{21}{x^2}\right) z \{ 0.1117 + z[4.421 + z(-9.207
        + 5.674\,z)]\} & \mathrm{for}\quad z > 0 \end{cases} \f] with \f[ z = \frac{x^2 -
        0.855}{x^2+3.42} \f]

        This approximation is accurate to within 3 per cent for all \f$x\f$ as long as
        \f$a<0.03\f$, which for the calculation of the Lyman-alpha cross section corresponds to a
        gas temperature above the cosmic microwave background temperature, i.e.
        \f$T_\mathrm{gas}>T_\mathrm{CMB}\f$. */
    double value(double a, double x);

    /** This function samples a random value from the probability distribution \f$P(u)\f$ defined
        by \f[ P(u) \propto \frac{\mathrm{e}^{-u^2}}{(u-x)^2+a^2} \f] where \f$a\f$ and \f$x\f$ are
        parameters and the proportionality factor is determined by normalization.

        We use the method described by Zheng et al. 2002 (ApJ, 578, 33-42, appendix), which is also
        used by Tasitsiomi 2006 (ApJ, 645, 792-813) and Laursen et al. 2009 (ApJ, 696, 853-869).
        The method is based on the rejection technique for sampling from a probability
        distribution. One could use the comparison function \f$g(u)\propto \left[ (u-x)^2 + a^2
        \right]^{-1}\f$ which can be integrated and inverted analytically so that it can be sampled
        using \f$u= x + a \tan \left[ \frac{\pi}{2}(2\mathcal{X}-1) \right] \f$, with
        \f$\mathcal{X}\f$ a uniform deviate.

        Because of the peculiar shape of \f$P(u)\f$, however, the comparison function \f$g(u)\f$ is
        better defined as \f[ g(u) \propto \begin{cases} \left[ (u-x)^2 + a^2 \right]^{-1} & u \le
        u_0 \\ \mathrm{e}^{-u_0^2}\left[ (u-x)^2 + a^2 \right]^{-1} & u > u_0 \end{cases} \f] where
        \f$u_0\f$ is determined by \f[ u_0 = \begin{cases} 0 & 0 \le x < 0.2 \\ x - 0.01 a^{1/6}
        \mathrm{e}^{1.2x} & 0.2 \le x < x_\mathrm{cw}(a) \\ 4.5 & x \ge x_\mathrm{cw}(a)
        \end{cases} \f] and the core/wing transition \f$x_\mathrm{cw}\f$ can be approximated by \f[
        x_\mathrm{cw}(a) = 1.59 - 0.6 \log_{10} a - 0.03 \log_{10}^2 a. \f]

        The acceptance fractions are required to be \f$\mathrm{e}^{-􏰜u^2}\f$ and
        \f$\mathrm{e}^{-􏰜u^2} / \mathrm{e}^{-􏰜u_0^2}\f$ in the regions \f$u\le 􏰠u_0\f$ and
        \f$u>u_0\f$, respectively. A first random number uniformly distributed between 0 and 1
        determines which region we use by comparing it with p, where \f[ p =
        \frac{\int_{-\infty}^{u_0} g(u)\,\mathrm{d}u} {\int_{-\infty}^{+\infty} g(u)\,\mathrm{d}u}
        = \left( \theta_0 + \frac{\pi}{2} \right) \left[ \left( 1-\mathrm{e}^{-u_0^2} \right)
        \theta_0 + \left( 1+\mathrm{e}^{-u_0^2} \right) \frac{\pi}{2} \right]^{-1}, \quad \theta_0
        = \arctan \frac{u_0 -x}{a} . \f] */
    double sample(double a, double x, std::function<double()> uniform);
}

////////////////////////////////////////////////////////////////////

#endif
