/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef SPECIALFUNCTIONS_HPP
#define SPECIALFUNCTIONS_HPP

#include "Basics.hpp"

////////////////////////////////////////////////////////////////////

/** This namespace contains a number of special mathematical functions. It includes the usual
    mathematical special functions such as the gamma function and elliptic integrals, but also a
    number of specific functions that often occur in the program. Many functions are copied or
    adapted from the famous Numerical Recipes books. The version used is the second edition of
    Numerical Recipes in C++ (Press et al. 2002). Other functions are taken from the web or are
    just implemented ad hoc. */
namespace SpecialFunctions
{
    /** This function returns the logarithm of the Gamma function, i.e. \f$\ln\Gamma(a)\f$. The
        implementation is taken from the Numerical Recipes book. */
    double lngamma(double a);

    /** This function returns the Gamma function \f$\Gamma(a)\f$. */
    double gamma(double a);

    /** This function returns the (regularized) incomplete gamma function \f[ P(a,x) =
        \frac{\gamma(a,x)}{\Gamma(a)} = \frac{1}{\Gamma(a)} \int_0^x {\text{e}}^{-t}\, t^{a-1}\,
        {\text{d}}t. \f] The implementation is adapted from the Numerical Recipes book. */
    double incompleteGamma(double a, double x);

    /** This function returns the regularized incomplete Beta function \f$I_x(a,b)\f$, defined as
        \f[ I_x(a,b) = \frac{B_x(a,b)}{B(a,b)} = \frac{1}{B(a,b)} \int_0^x
        t^{a-1}\,(1-t)^{b-1}\,dt, \f] with the parameters \f$a\f$ and \f$b\f$ positive real numbers
        and the argument \f$x\f$ between 0 and 1. The implementation uses the representation of
        this function in terms of continued fractions, and is adapted from the Numerical Recipes
        book. */
    double betaRegularized(double x, double a, double b);

    /** This function returns the complete Beta function \f$B(a,b)\f$, defined as \f[ B(a,b) =
        \int_0^1 t^{a-1}\,(1-t)^{b-1}\,dt, \f] with the two parameters \f$a\f$ and \f$b\f$ positive
        real numbers. The implementation is the trivial translation of the formula \f[ B(a,b) =
        \dfrac{\Gamma(a)\,\Gamma(b)}{\Gamma(a+b)} \f] */
    double beta(double a, double b);

    /** This function returns the incomplete Beta function \f$B_x(a,b)\f$, defined as \f[ B_x(a,b)
        = \int_0^x t^{a-1}\,(1-t)^{b-1}\,dt, \f] with the parameters \f$a\f$ and \f$b\f$ positive
        real numbers and the argument \f$x\f$ between 0 and 1. */
    double beta(double x, double a, double b);

    /** This function returns the modified Bessel function \f$I_1(x)\f$ of the first kind of
        integer order 1. The implementation is based on a polynomial approximation for small
        \f$x\f$ and the product of a polynomial approximation and the exponential factor at large
        \f$x\f$. The implementation is adapted from Numerical Recipes. */
    double bessi1(double x);

    /** This function returns the modified Bessel function \f$K_1(x)\f$ of the second kind of
        integer order 1. The implementation is based on a polynomial approximation for small
        \f$x\f$ and the product of a polynomial approximation and the exponential factor at large
        \f$x\f$. The implementation is adapted from Numerical Recipes. */
    double bessk1(double x);

    /** Dawson's integral, defined as \f[ D(x) = \exp(-x^2)\int_0^x \exp(t^2)\,dt.\f] The
        implementation is based on the series \f[ D(x) = \lim_{h\rightarrow0} \frac{1}{\sqrt{\pi}}
        \sum_{\text{$n$ odd}} \frac{e^{-(x-nh)^2}}{n}, \f] which guarantees an exponential increase
        in accuracy. The implementation is adapted from Numerical Recipes, where additional tricks
        are included to make the computation even faster. */
    double dawson(double x);

    /** This function returns Carlson's elliptic integral \f$R_F\f$, defined as \f[ R_F(x,y,z) =
        \frac{1}{2} \int_0^\infty\frac{dt}{\sqrt{(t+x)(t+y)(t+z)}}. \f] The three arguments
        \f$x\f$, \f$y\f$ and \f$z\f$ must be nonnegative and at most one of them is zero. The
        implementation of this function is based on the duplication theorem and is taken from
        Numerical Recipes. */
    double rf(double x, double y, double z);

    /** This function returns Carlson's elliptic integral \f$R_D\f$, defined as \f[ R_D(x,y,z) =
        \frac{3}{2} \int_0^\infty \frac{dt}{\sqrt{(t+x)(t+y)(t+z)^3}}. \f] The first two arguments
        \f$x\f$ and \f$y\f$ must be nonnegative and at most one of them is zero, whereas the last
        argument \f$z\f$ must be positive. The implementation, based on the duplication theorem, is
        taken from Numerical Recipes. */
    double rd(double x, double y, double z);

    /** This function returns Carlson's elliptic integral \f$R_J\f$, defined as \f[ R_J(x,y,z,p) =
        \frac{3}{2} \int_0^\infty \frac{dt}{(t+p)\sqrt{(t+x)(t+y)(t+z)}}. \f] The first three
        arguments \f$x\f$, \f$y\f$ and \f$z\f$ must be nonnegative and at most one of them is zero,
        whereas the last argument \f$p\f$ must be nonzero. The implementation, based on the
        duplication theorem, is taken from Numerical Recipes. */
    double rj(double x, double y, double z, double p);

    /** This function returns Carlson's degenerate elliptic integral \f$R_C\f$, defined as \f[
        R_C(x,y) = \frac{1}{2} \int_0^\infty \frac{dt}{(t+y)\sqrt{t+x}}. \f] Here \f$x\f$ must be
        nonnegative and \f$y\f$ must be nonzero. The implementation, based on the duplication
        theorem, is taken from Numerical Recipes. */
    double rc(double x, double y);

    /** This function returns the incomplete elliptic integral of the first kind, defined as \f[
        F(x,k) = \int_0^x \frac{dt}{\sqrt{(1-t^2)(1-k^2t^2)} }. \f] The arguments must satisfy
        \f$0\leq x\leq 1\f$ and \f$0\leq k \leq 1\f$. For the evaluation of this function, we
        express it in terms of Carlson's elliptic integral of the first kind, \f[ F(x,k) = x\,
        R_F(1-x^2,1-k^2x^2,1). \f] It is important to note that different incompatible argument
        conventions are being used for the elliptic integrals. Instead of the coordinate \f$x\f$
        one often uses the amplitude \f$\phi\f$, where both are related through \f$x=\sin\phi\f$.
        For the second parameter, one can use the modulus \f$k\f$ (used in SKIRT), the parameter
        \f$m\f$ or the modular angle \f$\alpha\f$. These are related through
        \f$m=k^2=\sin^2\alpha\f$. We use the conventions adopted in Maple. */
    double EllipticF(double x, double k);

    /** This function returns the complete elliptic integral of the first kind, defined as \f[ K(k)
        = \int_0^1 \frac{dt}{ \sqrt{(1-t^2)(1-k^2t^2)} }, \f] where \f$0\leq k\leq 1\f$. The
        function is no more than a function call to the incomplete elliptic integral of the first
        kind \f$F(x,k)\f$ with \f$x=1\f$. See also the note about argument conventions there. */
    double EllipticK(double k);

    /** This function returns the incomplete elliptic integral of the second kind, defined as \f[
        E(x,k) = \int_0^x \frac{\sqrt{1-k^2t^2}\, dt}{ \sqrt{1-t^2} }. \f] The arguments must
        satisfy \f$0\leq x\leq 1\f$ and \f$0\leq k \leq 1\f$. For the evaluation of this function,
        we express it in terms of Carlson's elliptic integrals of the first and second kind, \f[
        E(x,k) = x\, R_F(1-x^2,1-k^2x^2,1) - \frac{1}{3}\,x^3\,k^2\,R_D(1-x^2,1-k^2x^2,1). \f] It
        is important to note that different incompatible argument conventions are being used for
        the elliptic integrals. Instead of the coordinate \f$x\f$ one often uses the amplitude
        \f$\phi\f$, where both are related through \f$x=\sin\phi\f$. For the second parameter, one
        can use the modulus \f$k\f$ (used in SKIRT), the parameter \f$m\f$ or the modular angle
        \f$\alpha\f$. These are related through \f$m=k^2=\sin^2\alpha\f$. We use the conventions
        adopted in Maple. */
    double EllipticE(double x, double k);

    /** This function returns the complete elliptic integral of the second kind, defined as \f[
        E(k) = \int_0^1 \frac{\sqrt{1-k^2t^2}\, dt}{ \sqrt{1-t^2} },\f] where \f$0\leq k\leq 1\f$.
        The function is no more than a function call to the incomplete elliptic integral of the
        second kind \f$E(x,k)\f$ with \f$x=1\f$. See also the note about argument conventions
        there. */
    double EllipticE(double k);

    /** This function returns the incomplete elliptic integral of the third kind, defined as \f[
        \Pi(x,\nu,k) = \int_0^x \frac{dt}{(1-\nu t^2) \sqrt{(1-t^2) (1-k^2t^2)} }. \f] The
        arguments must satisfy \f$0\leq x\leq 1\f$ and \f$0\leq k \leq 1\f$. For the evaluation of
        this function, we express it in terms of Carlson's elliptic integrals of the first and
        third kind, \f[ \Pi(x,\nu,k) = x\, R_F(1-x^2,1-k^2x^2,1) +
        \frac{1}{3}\,\nu\,x^3\,R_J(1-x^2,1-k^2x^2,1,1-\nu\,x^2). \f] It is important to note that
        different incompatible argument conventions are being used for the elliptic integrals.
        Instead of the coordinate \f$x\f$ one often uses the amplitude \f$\phi\f$, where both are
        related through \f$x=\sin\phi\f$. For the second parameter, one can use the modulus \f$k\f$
        (used in SKIRT), the parameter \f$m\f$ or the modular angle \f$\alpha\f$. These are related
        through \f$m=k^2=\sin^2\alpha\f$. Finally, there are different sign conventions for the
        characteristic \f$\nu\f$, which is sometimes defined with a different signature (e.g.
        Abramowitz & Stegun 1965). We use the conventions adopted in Maple. */
    double EllipticPi(double x, double nu, double k);

    /** This function returns the complete elliptic integral of the third kind, defined as \f[
        \Pi(\nu,k) = \int_0^1 \frac{dt}{(1-\nu t^2) \sqrt{(1-t^2) (1-k^2t^2)} },\f] where \f$0\leq
        k\leq 1\f$. The function is no more than a function call to the incomplete elliptic
        integral of the third kind \f$\Pi(x,\nu,k)\f$ with \f$x=1\f$. See also the note about
        argument conventions there. */
    double EllipticPi(double nu, double k);

    /** This function returns the value of the function defined as \f[ X(s) = \begin{cases} \;
        \dfrac{1}{\sqrt{1-s^2}}\, \text{arccosh} \left(\dfrac{1}{s}\right) & \qquad\text{if
        $0<s<1$,} \\ \; \dfrac{1}{\sqrt{s^2-1}}\, \text{arccos} \left(\dfrac{1}{s}\right) &
        \qquad\text{if $s>1$.} \end{cases} \f] This function is completely smooth (it is the same
        function if one considers it in the complex plane). */
    double functionX(double s);

    /** This function is the Lambert \f$W\f$ function, also known as the product log function for
        the principle branch. The Lambert \f$W\f$ function is generally defined as the inverse of
        the function \f[ w \rightarrow f(w) = w\,e^w. \f] For real \f$z\f$ in the interval
        \f$[-1/e,0[\f$, the equation \f[ z=w\,e^w \f] has two real solutions, one smaller than and
        one greater than \f$-1\f$. These define the two only real branches of this complex
        function: the principle branch \f$W \equiv W_0\f$ returns a value \f$W(z)>-1\f$, whereas
        the other branch, labelled \f$W_{-1}\f$ returns a value \f$W_{-1}(z)<-1\f$. For more
        details on this function, see the help facility in Maple on the command <tt>LambertW</tt>
        or in Mathematica on the command <tt>ProductLog</tt>. The adopted implementation for the
        function \f$W_0(z)\f$ is a combination of the C-code available from the homepage of <a
        href="http://keithbriggs.info">Keith Briggs</a> and the code from the GNU scientific
        library <a href="https://www.gnu.org/software/gsl/">GSL</a>. */
    double LambertW(double z);

    /** This function is the Lambert \f$W\f$ function, also known as the product log function, for
        the branch \f$-1\f$. The Lambert \f$W\f$ function is generally defined as the inverse of
        the function \f[ w \rightarrow f(w) = w\,e^w. \f] For real \f$z\f$ in the interval
        \f$[-1/e,0[\f$, the equation \f[ z=w\,e^w \f] has two real solutions, one smaller than and
        one greater than \f$-1\f$. These define the two only real branches of this complex
        function: the principle branch \f$W \equiv W_0\f$ returns a value \f$W(z)>-1\f$, whereas
        the other branch, labelled \f$W_{-1}\f$ returns a value \f$W_{-1}(z)<-1\f$. For more
        details on this function, see the help facility in Maple on the command <tt>LambertW</tt>
        or in Mathematica on the command <tt>ProductLog</tt>. The adopted implementation for the
        function \f$W_{-1}(z)\f$ is a combination of the C-code available from the homepage of <a
        href="http://members.lycos.co.uk/keithmbriggs">Keith Briggs</a> and the code from the GNU
        scientific library <a href="http://sources.redhat.com/gsl">GSL</a>. The code returns an
        error message and the program is halted when the input value does not lie within the
        interval \f$[-1/e,0[\f$. Comparison of the results of this implementation with the results
        from the Maple function yields identical values. Only for values of \f$z\f$ extremely close
        to the limiting value zero, we get some minor discrepancies, e.g. the difference between
        the two results is 0.015% for \f$z=-10^{-9}\f$. */
    double LambertW1(double z);

    /** This function returns the Debye function of order \f$n\f$, defined as \f[ D_n(x) =
        \frac{n}{x^n} \int_0^x \frac{t^n\,{\text{d}}t}{{\text{e}}^t-1}. \f] The order \f$n\f$
        should be an integer number between 1 and 20; the argument \f$x\f$ is a positive real
        number. The implementation is adapted from Richard Mathar's website at MPIA Heidelberg. */
    double DebyeD(int n, double x);

    /** This function returns a generalized logarithmic function \f${\text{gln}}(p,x)\f$, defined
        for \f$x>0\f$ and arbitrary real \f$p\f$ as \f[ {\text{gln}}(p,x) = \int_1^x
        t^{-p}\,{\text{d}}t = \begin{cases} \; \dfrac{x^{1-p}-1}{1-p} & p\ne1 \\ \; \ln x & p=1
        \end{cases}. \f] This function is included in this library of special functions because the
        power law expression diverges as \f$p\f$ approaches 1. In that case we can use the
        expansion \f[ {\text{gln}}(p,x) \approx \ln x + \frac12\,(1-p)\ln^2x + \frac16\,(1-p)^2
        \ln^3x + \frac{1}{24}\,(1-p)^3\ln^4x + \ldots \f] */
    double gln(double p, double x);

    /** This function returns the difference between two values of the generalized logarithmic
        function \f${\text{gln}}(p,x)\f$ with the same exponent \f$p\f$. Compared to simply
        subtracting the two generalized logarithms, this function achieves much better accuracy for
        large arguments \f$x_1\f$ and \f$x_2\f$ by using the identity \f[ {\text{gln2}}(p,x_1,x_2)
        = {\text{gln}}(p,x_1) - {\text{gln}}(p,x_2) = (x_2)^{1-p} \,
        {\text{gln}}(p,\frac{x_1}{x_2}). \f] */
    double gln2(double p, double x1, double x2);

    /** This function returns a generalized exponential function \f${\text{gexp}}(p,x)\f$, defined
        as the inverse of the function gln. In formula it is defined as \f[ {\text{gexp}}(p,x) =
        \begin{cases} \; ((1-p)\,x+1)^{\frac{1}{1-p}} & p\ne1 \\ \; {\text{e}}^p & p=1 \end{cases}
        \f] When \f$p\f$ approaches 1, we use the expansion \f[ {\text{gexp}}(p,x) \approx
        {\text{e}}^x\left[ 1 - \frac12\,(1-p)\,x^2 + \frac{1}{24}\, (1-p)^2\, (3x+8)\,x^3 -
        \frac{1}{48}\, (1-p)^3\, (x^2+8x+12)\, x^4 + \ldots \right] \f] */
    double gexp(double p, double x);

    /** This function returns the logarithmic mean \f$M(x_1,x_2)\f$ of two nonnegative values
        \f$x_1\geq 0\f$ and \f$x_2\geq 0\f$, defined as \f[ M(x_1,x_2) = \begin{cases} 0 &
        \mathrm{if}\, x_1=0 \,\mathrm{or}\, x_2=0 \\ x_1 & \mathrm{if}\, x_1=x_2 \\ \dfrac{x_2 -
        x_1}{\ln x_2 - \ln x_1}&\mathrm{otherwise}.\end{cases} \f] The function is invariant for
        swapping of it arguments, and it can be proven that the logarithmic mean lies between the
        geometric mean and the arithmetic mean, \f[ \sqrt{(x_1 x_2)} \leq M(x_1,x_2) \leq
        \frac{x_1+x_2}{2}, \f] where the equality realizes for \f$x_1=x_2\f$.

        The function is implemented here because the quotient becomes numerically unstable when
        \f$x_1\approx x_2\f$. In that case, the function value is calculated using the substitution
        \f$x=\dfrac{x_2}{x_1}-1\f$ and the following expansion: \f[ M(x_1,x_2) = x_1 \,
        \frac{x}{\ln(1+x)} = \dfrac{x_1}{1 - \dfrac{x}{2} + \dfrac{x^2}{3} - \dfrac{x^3}{4} + ...}
        \f] */
    double lnmean(double x1, double x2);

    /** This function returns the logarithmic mean of two values as described for the two-argument
        lnmean() function in this class, given also the natural logarithm of these two values. In
        cases where these logarithms are available at the call site anyway, this function is more
        efficient than its two-argument equivalent. */
    double lnmean(double x1, double x2, double lnx1, double lnx2);
}

////////////////////////////////////////////////////////////////////

#endif
