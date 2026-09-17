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
        by \f[ P(u) \propto \frac{\mathrm{e}^{-u^2}}{(u-x)^2+a^2} \f] where \f$a\f$ and \f$x\f$ are
        parameters given as arguments and the proportionality factor is determined by
        normalization. The third argument specifies the random generator to be used by the
        function.

        We use the method described by Zheng et al. 2002 (ApJ, 578, 33-42, appendix), a variation
        of which is also used by many other authors including Tasitsiomi 2006 (ApJ, 645, 792-813),
        Laursen et al. 2009 (ApJ, 696, 853-869), Smith et al. 2015 (MNRAS, 449, 4336–4362) and
        Michel-Dansac et al. 2020 (A\&A). The method is based on the rejection technique for
        sampling from a probability distribution. One could use the comparison function
        \f$g(u)\propto \left[ (u-x)^2 + a^2 \right]^{-1}\f$ which can be integrated and inverted
        analytically so that it can be sampled using \f$u= x + a \tan \left[ \frac{\pi}{2}
        (2\mathcal{X}-1) \right] \f$, with \f$\mathcal{X}\f$ a uniform deviate.

        Because of the peculiar shape of \f$P(u)\f$, however, the comparison function \f$g(u)\f$ is
        better defined as \f[ g(u) \propto \begin{cases} \left[ (u-x)^2 + a^2 \right]^{-1} & u \le
        u_0 \\ \mathrm{e}^{-u_0^2}\left[ (u-x)^2 + a^2 \right]^{-1} & u > u_0 \end{cases} \f] where
        \f$u_0\f$ is determined following Michel-Dansac et al. 2020 (A\&A) by \f{align*} u_0 = &\;
        2.648963+2.014446 \zeta+0.351479\zeta^2 \\ &+ x (-4.058673-3.675859\zeta-0.640003\zeta^2 \\
        &+ x (3.017395+2.117133\zeta+0.370294\zeta^2 \\ &+ x
        (-0.869789-0.565886\zeta-0.096312\zeta^2 \\ &+ x (0.110987+0.070103\zeta+0.011557\zeta^2 \\
        &+ x (-0.005200-0.003240\zeta-0.000519\zeta^2))))) \f} where \f$\zeta = \log_{10}(a)\f$.

        The acceptance fractions are required to be \f$\mathrm{e}^{-u^2}\f$ and
        \f$\mathrm{e}^{-u^2} / \mathrm{e}^{-u_0^2}\f$ in the regions \f$u\le u_0\f$ and
        \f$u>u_0\f$, respectively. A first random number uniformly distributed between 0 and 1
        determines which region we use by comparing it with p, where \f[ p =
        \frac{\int_{-\infty}^{u_0} g(u)\,\mathrm{d}u} {\int_{-\infty}^{+\infty} g(u)\,\mathrm{d}u}
        = \left( \theta_0 + \frac{\pi}{2} \right) \left[ \left( 1-\mathrm{e}^{-u_0^2} \right)
        \theta_0 + \left( 1+\mathrm{e}^{-u_0^2} \right) \frac{\pi}{2} \right]^{-1}, \quad \theta_0
        = \arctan \frac{u_0 -x}{a} . \f]

        Finally, for larger values of \f$x\f$, the distribution \f$P(u)\f$ can successfully be
        approximated by a Gaussian distribution centered on \f$1/x\f$. Following Smith et al. 2015
        (MNRAS, 449, 4336–4362) and Michel-Dansac et al. 2020 (A\&A), we use this approximation
        for \f$x \ge 8\f$. */
    double sample(double a, double x, Random* random);
}

////////////////////////////////////////////////////////////////////

#endif
