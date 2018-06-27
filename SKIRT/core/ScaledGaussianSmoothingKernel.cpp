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
    const double N = 2.56810060330949540082;    // front factor N
    const double A = -5.85836755024609305208;   // exponent factor -1/(2 sigma^2)

    // the number of equidistant points in the cumulative distribution grid
    const int Nu = 400;
}

//////////////////////////////////////////////////////////////////////

void ScaledGaussianSmoothingKernel::setupSelfBefore()
{
    SmoothingKernel::setupSelfBefore();

    double du = 1.0/Nu;
    _Xv.resize(Nu+1);
    _Xv[0] = 0.;
    for (int k=1; k<Nu; k++)
    {
        double u = k*du;
        _Xv[k] = _Xv[k-1] + ( 4.*M_PI * density(u) * u*u * du );
    }
    _Xv[Nu] = 1.;
}

//////////////////////////////////////////////////////////////////////

double ScaledGaussianSmoothingKernel::density(double u) const
{
    if (u < 0. || u > 1.) return 0.;
    return N * exp(A*u*u);
}

//////////////////////////////////////////////////////////////////////

double ScaledGaussianSmoothingKernel::generateRadius() const
{
    double X = random()->uniform();
    int k = NR::locateClip(_Xv,X);
    double p = (X-_Xv[k])/(_Xv[k+1]-_Xv[k]);
    double u = (k+p)/Nu;
    return u;
}

//////////////////////////////////////////////////////////////////////
