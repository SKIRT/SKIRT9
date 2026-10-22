/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "VoigtProfile.hpp"
#include "FatalError.hpp"
#include "Random.hpp"
#include <complex>

////////////////////////////////////////////////////////////////////

double VoigtProfile::value(double a, double x)
{
    // Humlíček (1982) approximation for w(z = x + i*a), valid for all x and a >= 0
    std::complex<double> t(a, -x);

    double s = std::fabs(x) + a;

    if (s >= 15.0)
    {
        // Region I
        return std::real(t * 0.5641896 / (0.5 + t * t));
    }
    else if (s >= 5.5)
    {
        // Region II
        std::complex<double> u = t * t;
        return std::real(t * (1.410474 + u * 0.5641896) / (0.75 + u * (3.0 + u)));
    }
    else if (a >= 0.195 * std::fabs(x) - 0.176)
    {
        // Region III
        return std::real((16.4955 + t * (20.20933 + t * (11.96482 + t * (3.778987 + t * 0.5642236))))
                         / (16.4955 + t * (38.82363 + t * (39.27121 + t * (21.69274 + t * (6.699398 + t))))));
    }
    else
    {
        // Region IV
        std::complex<double> u = t * t;
        return std::real(
            std::exp(u)
            - t
                  * (36183.31
                     - u
                           * (3321.9905
                              - u * (1540.787 - u * (219.0313 - u * (35.76683 - u * (1.320522 - u * 0.56419))))))
                  / (32066.6
                     - u
                           * (24322.84
                              - u
                                    * (9022.228
                                       - u * (2186.181 - u * (364.2191 - u * (61.57037 - u * (1.841439 - u))))))));
    }
}

////////////////////////////////////////////////////////////////////

namespace
{
    // the maximum number of attempts in the rejection method loops
    constexpr int maxAttempts = 100000;

    // sample from the distribution for x >= 0 using a piecewise Lorentzian comparison function
    double sampleLorentzian(double a, double x, Random* random)
    {
        // determine the separation u0 that (approximately) maximizes the acceptance rate; any u0 >= 0
        // yields an exact sampling procedure, the choice merely affects the efficiency; minimizing the
        // normalization of the comparison function shows that the distance d = x - u0 approximately
        // satisfies d^2 exp(2xd) = a exp(x^2) / (2 pi x), so that xd = W(y) with W the Lambert function
        // and y = x sqrt(a exp(x^2) / (2 pi x)); because high accuracy is not needed, we use the
        // approximation W(y) = L (1 - ln(1+L) / (2+L)) with L = ln(1+y) (Winitzki 2003), evaluated
        // from ln(y) to avoid overflow; the result is kept strictly below x because for extremely small
        // a the distance may be smaller than the spacing between floating point numbers near x
        double u0 = 0.;
        if (x > 0.)
        {
            double lny = 0.5 * (log(a * x / (2. * M_PI)) + x * x);
            double L = lny > 0. ? lny + log1p(exp(-lny)) : log1p(exp(lny));
            double W = L * (1. - log1p(L) / (2. + L));
            u0 = max(0., min(x - W / x, nextafter(x, 0.)));
        }

        // the weights of the Lorentzian to the left and right of u0, and the probability for the left side
        double wl = atan2(a, x - u0);
        double wr = atan2(a, u0 - x);
        double p = wl / (wl + exp(-u0 * u0) * wr);

        // perform the rejection method loop for a maximum number of attempts
        int n = maxAttempts;
        while (n--)
        {
            // generate a random sample from the selected part of the comparison function,
            // inverting the Lorentzian cumulative distribution in a numerically stable form
            double u;
            if (random->uniform() <= p)
                u = x - a / tan(wl * random->uniform());
            else
                u = x + a / tan(wr * random->uniform());

            // determine the acceptance fraction, which never exceeds unity because u0 >= 0
            double fraction = u > u0 ? exp((u0 - u) * (u0 + u)) : exp(-u * u);

            // accept or reject the sample
            if (random->uniform() < fraction) return u;
        }

        // if none of the attempts were accepted, abort
        throw FATALERROR("Sampling from Voigt profile has failed");
    }

    // sample from the distribution for x >= 0 and rho = sqrt(x^2+a^2) >= 8 using a Gaussian
    // comparison function; the caller must verify that the narrow peak at u = x (if any) is negligible
    double sampleGaussian(double a, double x, double rho, Random* random)
    {
        // the Gaussian proposal matching the distribution to second order in u, expressed in
        // scaled quantities so that nothing overflows for extremely large x or a
        double xs = x / rho;
        double as = a / rho;
        double alpha = 1. + (as - xs) * (as + xs) / (rho * rho);
        double beta = 2. * xs / rho;
        double mean = 0.5 * beta / alpha;
        double sigma = 1. / sqrt(2. * alpha);

        // the ratio of the target distribution and the proposal, normalized to unity at u = 0
        auto ratio = [=](double u) {
            double us = u / rho - xs;
            return exp((alpha - 1.) * u * u - beta * u) / (us * us + as * as);
        };

        // an upper bound for the ratio over the range |u| <= 5.5, outside of which the proposal has a
        // negligible weight (below 1e-14); the maximum is attained at or very close to the edges of
        // this range, and the factor 1.1 provides a safety margin (verified numerically for all a, x)
        double bound = 1.1 * max(ratio(-5.5), ratio(5.5));

        // perform the rejection method loop for a maximum number of attempts
        int n = maxAttempts;
        while (n--)
        {
            double u = mean + sigma * random->gauss();
            if (random->uniform() * bound < ratio(u)) return u;
        }

        // if none of the attempts were accepted, abort
        throw FATALERROR("Sampling from Voigt profile has failed");
    }
}

////////////////////////////////////////////////////////////////////

double VoigtProfile::sample(double a, double x, Random* random)
{
    // the distribution is symmetric in x, so sample for |x| and restore the sign
    double ax = abs(x);

    // far from the line core, use the Gaussian comparison function unless the narrow peak at u = x,
    // which it ignores, carries a non-negligible weight (relative weight below 1e-15); this
    // exception only occurs for extremely small values of a (below about 1e-10)
    double rho = hypot(ax, a);
    bool gaussian = rho >= 8. && (a >= ax || ax * ax - log(rho * rho / a) > 35.);

    double u = gaussian ? sampleGaussian(a, ax, rho, random) : sampleLorentzian(a, ax, random);
    return x < 0. ? -u : u;
}

////////////////////////////////////////////////////////////////////
