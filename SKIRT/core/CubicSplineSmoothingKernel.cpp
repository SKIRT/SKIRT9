/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "CubicSplineSmoothingKernel.hpp"
#include "NR.hpp"
#include "Random.hpp"

//////////////////////////////////////////////////////////////////////

namespace
{
    const int _Nu = 400;
}

//////////////////////////////////////////////////////////////////////

void CubicSplineSmoothingKernel::setupSelfBefore()
{
    SmoothingKernel::setupSelfBefore();

    double du = 1.0 / _Nu;
    _Xv.resize(_Nu + 1);
    for (int k = 0; k <= _Nu; k++)
    {
        double u = k * du;
        double u2 = u * u;
        double u3 = u * u2;
        if (u < 0.5)
            _Xv[k] = u3 * (32.0 / 3.0 - 192.0 / 5.0 * u2 + 32.0 * u3);
        else
            _Xv[k] = -1.0 / 15.0 - 64.0 * u3 * (-1.0 / 3.0 + 0.75 * u - 0.6 * u2 + u3 / 6.0);
    }
}

//////////////////////////////////////////////////////////////////////

double CubicSplineSmoothingKernel::density(double u) const
{
    if (u < 0.0 || u >= 1.0)
        return 0.0;
    else if (u < 0.5)
        return 8.0 / M_PI * (1.0 - 6.0 * u * u * (1.0 - u));
    else
        return 8.0 / M_PI * 2.0 * (1.0 - u) * (1.0 - u) * (1.0 - u);
}

//////////////////////////////////////////////////////////////////////

double CubicSplineSmoothingKernel::columnDensity(double q) const
{
    if (q < 0.0 || q >= 1.0) return 0.0;

    double q2 = q * q;
    if (q < 1e-6) return 6.0 / M_PI * (1.0 - 8.0 * M_LN2 * q2);

    double s = sqrt((1.0 - q) * (1.0 + q));
    if (q < 0.5)
    {
        double t = sqrt((1.0 - 2.0 * q) * (1.0 + 2.0 * q));
        double p1 = (4.0 + 26.0 * q2) * s;
        double p2 = (1.0 + 26.0 * q2) * t;
        double p3 = 18.0 * q2 * q2 * log(2.0 * q / (1.0 + t));
        double p4 = 6.0 * q2 * (4.0 + q2) * log(2.0 * (1.0 + s) / (1.0 + t));
        return 2.0 / M_PI * (p1 - p2 - p3 - p4);
    }
    else
    {
        double p1 = (2.0 + 13.0 * q2) * s;
        double p2 = 3.0 * q2 * (4.0 + q2) * log(q / (1.0 + s));
        return 4.0 / M_PI * (p1 + p2);
    }
}

//////////////////////////////////////////////////////////////////////

double CubicSplineSmoothingKernel::generateRadius() const
{
    double X = random()->uniform();
    int k = NR::locateClip(_Xv, X);
    double p = (X - _Xv[k]) / (_Xv[k + 1] - _Xv[k]);
    double u = (k + p) / _Nu;
    return u;
}

//////////////////////////////////////////////////////////////////////
