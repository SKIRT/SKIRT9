/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "VoigtProfile.hpp"
#include "FatalError.hpp"

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

double VoigtProfile::sample(double a, double x, std::function<double()> uniform)
{
    // make x positive and remember the orginal sign
    double sign = 1.;
    if (x < 0.)
    {
        sign = -1.;
        x = -x;
    }

    // determine the core/wing transition corresponding to a
    double loga = log10(a);
    double xcw = 1.59 - 0.60 * loga - 0.03 * loga * loga;

    // determine the comparison function separation corresponding to a and x
    double u0 = 0.;
    if (x >= xcw)
        u0 = 4.5;
    else if (x >= 0.2)
        u0 = x - 0.01 * pow(a, 1. / 6.) * exp(1.2 * x);

    // calculate the cumulative separation point
    double theta0 = atan((u0 - x) / a);
    double p = (theta0 + M_PI_2) / ((1. - exp(-u0 * u0)) * theta0 + (1. + exp(-u0 * u0)) * M_PI_2);

    // perform the rejection method loop for a finite number of times
    int n = 1000000;
    while (n--)
    {
        // determine which one of the two comparison functions to use
        double left, right;
        if (uniform() <= p)
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
        double u = x + a * tan((right - left) * uniform() + left);

        // determine the acceptance/rejection fraction
        double fraction = exp(-u * u);
        if (u > u0) fraction /= exp(-u0 * u0);

        // accept or reject the sample
        if (uniform() < fraction) return u * sign;
    }

    // if the loop fails, abort
    throw FATALERROR("Sampling from Voigt profile has failed");
}

////////////////////////////////////////////////////////////////////
