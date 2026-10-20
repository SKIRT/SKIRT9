/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef VOIGTPROFILE_HPP
#define VOIGTPROFILE_HPP

#include "Basics.hpp"
class Random;

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
        parameter and the function is normalized so that \f$H(a,0) \approx 1\f$ for \f$a\ll 1\f$.
        For reference, the following table lists \f$H(a,0)\f$ for a few values of \f$a\f$.

        | \f$a\f$ | \f$H(a,0)\f$ |
        |---------|--------------|
        | 0.03    | 0.967029     |
        | 0.01    | 0.988815     |
        | 0.001   | 0.998873     |
        | 0.0001  | 0.999887     |

        We use the rational-approximation method of Humlíček 1982 (JQSRT, 27, 437-444) for the
        related complex probability function (Faddeeva function) \f$w(z) = \mathrm{e}^{-z^2}
        \mathrm{erfc}(-iz)\f$, which relates to \f$H(a,x)\f$ through \f$H(a,x) =
        \mathrm{Re}\,w(x+ia)\f$. The method evaluates \f$w(z)\f$ using one of four rational
        approximations, selected depending on the magnitude of \f$x\f$ and \f$a\f$: an inexpensive
        asymptotic form when \f$z=x+ia\f$ lies far from the origin, and successively higher-order
        forms as it approaches the origin, with a dedicated approximation for the region close to
        the real axis (small \f$a\f$) where a single rational approximation is known to lose
        accuracy. Overall the method is accurate to about four significant digits, substantially
        better than the approximation of Smith et al. 2015 (MNRAS, 449, 4336–4362) used
        previously, whose error can reach 1 per cent for \f$a<0.03\f$. */
    double value(double a, double x);

    /** This function samples a random value from the probability distribution \f$P(u)\f$ defined
        by \f[ P(u) \propto \frac{\mathrm{e}^{-u^2}}{(u-x)^2+a^2} \f] where \f$a>0\f$ and \f$x\f$ are
        parameters given as arguments and the proportionality factor is determined by
        normalization. The third argument specifies the random generator to be used by the
        function. Because the distribution is symmetric under \f$(x,u)\to(-x,-u)\f$, we sample for
        \f$|x|\f$ and restore the sign afterwards. The function uses one of two exact rejection
        methods, depending on the values of \f$a\f$ and \f$x\f$.

        <b>Piecewise Lorentzian comparison function</b>

        The first method, used close to the line core, follows Zheng et al. 2002 (ApJ, 578, 33-42,
        appendix), a variation of which is also used by many other authors including Tasitsiomi 2006
        (ApJ, 645, 792-813), Laursen et al. 2009 (ApJ, 696, 853-869), Smith et al. 2015 (MNRAS, 449,
        4336–4362) and Michel-Dansac et al. 2020 (A\&A). The comparison function is \f[ g(u) \propto
        \begin{cases} \left[ (u-x)^2 + a^2 \right]^{-1} & u \le u_0 \\ \mathrm{e}^{-u_0^2}\left[
        (u-x)^2 + a^2 \right]^{-1} & u > u_0 \end{cases} \f] which can be sampled by inverting the
        Lorentzian cumulative distribution on either side of \f$u_0\f$. The acceptance fractions
        are \f$\mathrm{e}^{-u^2}\f$ and \f$\mathrm{e}^{u_0^2-u^2}\f$ for \f$u\le u_0\f$ and
        \f$u>u_0\f$, respectively. The second fraction exceeds unity for \f$|u|<|u_0|\f$ when
        \f$u_0<0\f$, so the method is exact only if \f$u_0 \ge 0\f$. With the Lorentzian weights
        \f$w_\mathrm{l} = \arctan\left[a/(x-u_0)\right]\f$ and \f$w_\mathrm{r} = \pi -
        w_\mathrm{l}\f$ on either side of \f$u_0\f$ (evaluated with \f$\mathrm{atan2}\f$ to remain
        accurate for very small \f$a\f$), the left side is selected with probability \f$p =
        w_\mathrm{l} / C(u_0)\f$, where \f$C(u_0) = w_\mathrm{l} + \mathrm{e}^{-u_0^2}
        w_\mathrm{r}\f$ is the normalization of the comparison function.

        Any \f$u_0 \ge 0\f$ yields an exact procedure, but the acceptance rate is inversely
        proportional to \f$C(u_0)\f$. Previous versions of this function used a fitting formula for
        \f$u_0\f$ (Michel-Dansac et al. 2020) that was calibrated for neutral hydrogen. That formula
        produces \f$u_0<0\f$ for \f$3\times10^{-4} \lesssim a \lesssim 6\times10^{-3}\f$ near
        \f$x=0\f$, causing a bias of up to five per cent in the sampled distribution, and it
        produces acceptance rates below \f$10^{-7}\f$ for \f$a \lesssim 10^{-6}\f$ (as occurs, for
        example, for forbidden or intercombination lines of ions), causing the sampling procedure to
        fail. Instead, we now use a value of \f$u_0\ge0\f$ that approximately minimizes
        \f$C(u_0)\f$. Setting the derivative of \f$C(u_0)\f$ to zero and neglecting terms that are
        small for \f$a \ll 1\f$ shows that the distance \f$d=x-u_0\f$ satisfies \f$d^2\,
        \mathrm{e}^{2xd} = a\,\mathrm{e}^{x^2}/(2\pi x)\f$, with solution \f$d = W(y)/x\f$ where
        \f$W\f$ is the Lambert function and \f$y = x\sqrt{a\,\mathrm{e}^{x^2}/(2\pi x)}\f$.
        Because high accuracy is not needed, we use the approximation \f$W(y) \approx L\,[1 -
        \ln(1+L)/(2+L)]\f$ with \f$L=\ln(1+y)\f$ (Winitzki 2003, LNCS 2667, 780-789). The
        resulting value \f$u_0 = x - d\f$ is limited to the interval \f$[0, x[\f$. For most
        parameter values, the acceptance rate is within a few per cent of the rate obtained with the
        exact optimum, and it is never below about 55 per cent of that rate. The mean number of
        attempts per sample is below 25 for \f$a \ge 10^{-6}\f$, below 90 for \f$a \ge
        10^{-20}\f$, and below 150 for \f$a \ge 10^{-30}\f$.

        <b>Gaussian comparison function</b>

        Far from the line core, i.e. for \f$\rho = \sqrt{x^2+a^2} \ge 8\f$, the distribution is
        close to a Gaussian. Expanding the logarithm of the Lorentzian factor to second order in
        \f$u\f$ yields the comparison function \f[ q(u) \propto \exp\left( -\alpha u^2 + \beta u
        \right), \quad \alpha = 1 + \frac{a^2-x^2}{\rho^4}, \quad \beta = \frac{2x}{\rho^2}, \f] a
        Gaussian with mean \f$\beta/2\alpha\f$ (approximately \f$1/x\f$ for \f$a \ll x\f$) and
        variance \f$1/2\alpha\f$. A sample drawn from \f$q(u)\f$ is accepted with probability
        \f$r(u)/M\f$, where \f$r(u) = P(u)/q(u)\f$ normalized to \f$r(0)=1\f$, and \f$M = 1.1
        \max[r(-5.5), r(5.5)]\f$ is an upper bound for \f$r(u)\f$ over the range \f$|u| \le 5.5\f$
        (verified numerically). Outside of this range the comparison function has a negligible
        weight (below \f$10^{-14}\f$), so that the method is exact for all practical purposes. The
        mean number of attempts per sample is at most 1.8. This method ignores the narrow peak of
        \f$P(u)\f$ at \f$u=x\f$, which has relative weight \f$\approx \sqrt{\pi}\,\rho^2
        \mathrm{e}^{-x^2}/a\f$. Therefore it is used only if \f$a\ge x\f$ or if this weight is
        below \f$10^{-15}\f$; otherwise the first method is used.

        Previous versions of this function approximated the distribution for \f$x \ge 8\f$ by a
        Gaussian with mean \f$1/x\f$ and variance \f$1/2\f$. That approximation ignores the
        dependence on \f$a\f$; its cumulative distribution differs from the actual one by about
        0.3 per cent for small \f$a\f$ and by several per cent for \f$a \gtrsim 10\f$ (as occurs,
        for example, for some ions at low temperatures). */
    double sample(double a, double x, Random* random);
}

////////////////////////////////////////////////////////////////////

#endif
