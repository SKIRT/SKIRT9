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

double VoigtProfile::sample(double a, double x, Random* random)
{
    // make x positive and remember the orginal sign
    double sign = 1.;
    if (x < 0.)
    {
        sign = -1.;
        x = -x;
    }

    // when x is large, the distribution is essentially a Gaussian centered on 1/x
    if (x >= 8.) return sign / x + M_SQRT1_2 * random->gauss();

    // determine the comparison function separation corresponding to a and x
    double z = log10(a);
    double z2 = z * z;
    double u0 = 2.648963 + 2.014446 * z + 0.351479 * z2
                + x
                      * (-4.058673 - 3.675859 * z - 0.640003 * z2
                         + x
                               * (3.017395 + 2.117133 * z + 0.370294 * z2
                                  + x
                                        * (-0.869789 - 0.565886 * z - 0.096312 * z2
                                           + x
                                                 * (0.110987 + 0.070103 * z + 0.011557 * z2
                                                    + x * (-0.005200 - 0.003240 * z - 0.000519 * z2)))));

    // calculate the cumulative separation point
    double theta0 = atan((u0 - x) / a);
    double p = (theta0 + M_PI_2) / ((1. - exp(-u0 * u0)) * theta0 + (1. + exp(-u0 * u0)) * M_PI_2);

    // perform the rejection method loop for a maximum number of attempts
    int n = 100000;  // not the most efficient for certain a,x combinations (eg. Zn+29 at T=3K)
    while (n--)
    {
        // determine which one of the two comparison functions to use
        double left, right;
        if (random->uniform() <= p)
        {
            left = -M_PI_2;
            right = theta0;
        }
        else
        {
            left = theta0;
            right = M_PI_2;
        }

        // generate a random sample from the selected comparison function
        double u = x + a * tan((right - left) * random->uniform() + left);

        // determine the acceptance/rejection fraction
        double fraction = exp(-u * u);
        if (u > u0) fraction /= exp(-u0 * u0);

        // accept or reject the sample
        if (random->uniform() < fraction) return u * sign;
    }

    // if none of the attempts were accepted, abort
    throw FATALERROR("Sampling from Voigt profile has failed");
}

////////////////////////////////////////////////////////////////////
