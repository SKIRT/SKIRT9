/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef SERSICFUNCTION_HPP
#define SERSICFUNCTION_HPP

#include "Array.hpp"

//////////////////////////////////////////////////////////////////////

/** The SersicFunction class represents the Sersic function \f${\cal{S}}_n(s)\f$ of Sersic index
    \f$n\f$; see Sersic (1963). This function represents the dimensionless three-dimensional
    density of a normalized Sersic surface brightness distribution, \f[ I_n(s_{\text{p}}) =
    \frac{b_n^{2n}}{\pi\,\Gamma(2n+1)}\, {\text{e}}^{-b_ns_{\text{p}}^{1/n}}, \f] where \f$b_n\f$
    is a dimensionless constant that is defined by the condition that \f[ \pi \int_0^1
    I_n(s_{\text{p}})\, s_{\text{p}}\, {\text{d}}s_{\text{p}} = \frac12.\f] Its value can
    conveniently be approximated as \f[ b_n = 2n -\frac{1}{3} + \frac{4}{405n} +
    \frac{46}{25515n^2} + \frac{131}{1148175n^3}, \f] as derived by Ciotti & Bertin (1999, A&A 352,
    447–451). Given \f$I_n(s_{\text{p}})\f$, the Sersic function \f${\cal{S}}_n(s)\f$ can be
    determined by solving the Abel integral equation \f[ 2 \int_{s_{\text{p}}}^\infty \frac{
    {\cal{S}}_n(s)\, s\, {\text{d}}s }{ \sqrt{s^2-s_{\text{p}}^2} } = I_n(s_{\text{p}}). \f] The
    solution of this integral equation can be written as \f[ {\cal{S}}_n(s) = -\frac{1}{\pi}
    \int_s^\infty \frac{ \text{d}I_n }{ {\text{d}} s_{\text{p}}} (s_{\text{p}})\, \frac{ {\text{d}
    s_{\text{p}}}} {\sqrt{s_{\text{p}}^2-s^2}}. \f] As a C++ class, the SersicFunction class
    contains a vector with the values of \f${\cal{S}}_n(s)\f$ on a grid, distributed
    logarithmically between \f$s=10^{-6}\f$ and \f$s=10^4\f$. Also the cumulative mass function,
    \f[ {\cal{M}}_n(s) = 4\pi \int_0^s {\cal{S}}_n(s)\,s^2\, {\text{d}}s \f] is stored on the same
    grid. */
class SersicFunction
{
public:
    /** Constructor for the SersicFunction class. It creates a logarithmic grid of
        \f$s\f$-points ranging from \f$10^{-6}\f$ to \f$10^4\f$ and computes the integral for
        \f${\cal{S}}_n(s)\f$ on each of the grid points. An additional integration (with
        logarithmic interpolation between the grid points) is used to calculate the cumulative mass
        function. The accuracy of the integrations has been checked for the \f$n=4\f$ case through
        comparison with the tabulated values of Young (1976). For the case \f$n=1\f$ the
        integrations are checked directly against the analytical value: in this case we find \f[
        {\cal{S}}_1(s) = \frac{b_1^3}{2\pi^2}\, K_0(b_1s), \f] with \f$b_1 = 1.67834699\f$ and
        \f$K_0(x)\f$ the modified Bessel function of the second kind. In both cases, a relative
        accuracy better than 0.01% is typically reached. */
    SersicFunction(double n);

    /** The function returns the Sersic function \f${\cal{S}}_n(s)\f$ for the dimensionless radius
        \f$s\f$. Its value is determined by logarithmic interpolation on the internally stored
        grid. */
    double operator()(const double s) const;

    /** The function returns the cumulative mass function \f[ {\cal{M}}_n(s) = 4\pi \int_0^s
        {\cal{S}}_n(s')\, s'^2\, {\text{d}}s' \f] at the dimensionless radius \f$s\f$. Its value is
        determined by logarithmic interpolation on the internally stored grid. */
    double mass(const double s) const;

    /** The function returns the inverse function of the cumulative mass function; it solves the
        equation \f[ 4\pi \int_0^s {\cal{S}}_n(s')\, s'^2\, {\text{d}}s' = M \f] for \f$s\f$, given
        the normalized cumulative mass \f$M\f$, assumed to lie between 0 and 1. The value of
        \f$s\f$ is determined by logarithmic interpolation on the internally stored grid. */
    double inverseMass(const double M) const;

private:
    Array _sv;
    Array _Sv;
    Array _Mv;
};

//////////////////////////////////////////////////////////////////////

#endif
