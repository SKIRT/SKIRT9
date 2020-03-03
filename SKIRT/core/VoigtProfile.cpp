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
    // coefficients for the approximation function (Table A1, Smith+15)
    constexpr double A0 = 15.75328153963877;
    constexpr double A1 = 286.9341762324778;
    constexpr double A2 = 19.05706700907019;
    constexpr double A3 = 28.22644017233441;
    constexpr double A4 = 9.526399802414186;
    constexpr double A5 = 35.29217026286130;
    constexpr double A6 = 0.8681020834678775;
    constexpr double B0 = 0.0003300469163682737;
    constexpr double B1 = 0.5403095364583999;
    constexpr double B2 = 2.676724102580895;
    constexpr double B3 = 12.82026082606220;
    constexpr double B4 = 3.21166435627278;
    constexpr double B5 = 32.032981933420;
    constexpr double B6 = 9.0328158696;
    constexpr double B7 = 23.7489999060;
    constexpr double B8 = 1.82106170570;

    // calculation of the approximation (Appendix A1, Smith+15)
    double z = x * x;
    if (z <= 3.0) return exp(-z) * (1.0 - a * (A0 + A1 / (z - A2 + A3 / (z - A4 + A5 / (z - A6)))));
    if (z < 25.0) return exp(-z) + a * (B0 + B1 / (z - B2 + B3 / (z + B4 + B5 / (z - B6 + B7 / (z - B8)))));
    return 0.5 * M_2_SQRTPI * a / (z - 1.5 - 1.5 / (z - 3.5 - 5.0 / (z - 5.5)));
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
