/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "ScaledGaussianSmoothingKernel.hpp"
#include "NR.hpp"
#include "Random.hpp"

//////////////////////////////////////////////////////////////////////

namespace
{
    // the constants in the density formula calculated to high accuracy
    constexpr double N = 2.56810060330949540082;      // front factor N
    constexpr double sigma = 0.29214381374061638716;  // dispersion sigma
    constexpr double A = -5.85836755024609305208;     // exponent factor -1/(2 sigma^2)
    constexpr double B = 1.88060965502447027232;      // front factor column density N sqrt(2 pi) sigma

    // the number of equidistant points in the cumulative distribution grid
    const int Nu = 400;
}

//////////////////////////////////////////////////////////////////////

void ScaledGaussianSmoothingKernel::setupSelfBefore()
{
    SmoothingKernel::setupSelfBefore();

    double du = 1.0 / Nu;
    _Xv.resize(Nu + 1);
    _Xv[0] = 0.;
    for (int k = 1; k < Nu; k++)
    {
        double u = k * du;
        _Xv[k] = _Xv[k - 1] + (4. * M_PI * density(u) * u * u * du);
    }
    _Xv[Nu] = 1.;
}

//////////////////////////////////////////////////////////////////////

double ScaledGaussianSmoothingKernel::density(double u) const
{
    if (u < 0. || u > 1.) return 0.;
    return N * exp(A * u * u);
}

//////////////////////////////////////////////////////////////////////

double ScaledGaussianSmoothingKernel::columnDensity(double q) const
{
    if (q < 0. || q > 1.) return 0.;
    double s = sqrt((1.0 - q) * (1.0 + q));
    return B * exp(A * q * q) * erf(s / (M_SQRT2 * sigma));
}

//////////////////////////////////////////////////////////////////////

double ScaledGaussianSmoothingKernel::generateRadius() const
{
    double X = random()->uniform();
    int k = NR::locateClip(_Xv, X);
    double p = (X - _Xv[k]) / (_Xv[k + 1] - _Xv[k]);
    double u = (k + p) / Nu;
    return u;
}

//////////////////////////////////////////////////////////////////////
