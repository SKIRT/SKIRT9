/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "VoigtProfile.hpp"
#include "FatalError.hpp"
#include "Random.hpp"

////////////////////////////////////////////////////////////////////

double VoigtProfile::value(double a, double x)
{
    double x2 = x * x;
    double H = exp(-x2);
    double z = (x2 - 0.855) / (x2 + 3.42);
    if (z > 0.)
    {
        H += 0.5 * M_2_SQRTPI * a / (x2 + 1.) * (1. + 21. / x2) * z * (0.1117 + z * (4.421 + z * (-9.207 + 5.674 * z)));
    }
    return H;
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
    int n = 10000;
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
