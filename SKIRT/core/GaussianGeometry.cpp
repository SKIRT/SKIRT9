/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "GaussianGeometry.hpp"
#include "Random.hpp"

//////////////////////////////////////////////////////////////////////

void GaussianGeometry::setupSelfBefore()
{
    SpheGeometry::setupSelfBefore();

    // calculate cached values
    _rho0 = 1.0 / pow(sqrt(2.0 * M_PI) * _sigma, 3);

    // grid with values of the cumulative mass
    int N = 401;
    double logtmin = -4.0;
    double logtmax = 4.0;
    double dlogt = (logtmax - logtmin) / (N - 1.0);
    _rv.resize(N);
    _Xv.resize(N);
    _Xv[0] = 0.0;
    for (int i = 1; i < N - 1; i++)
    {
        double logt = logtmin + i * dlogt;
        double t = pow(10.0, logt);
        _rv[i] = M_SQRT2 * _sigma * t;  // t = r / sigma / sqrt(2)
        _Xv[i] = erf(t) - M_2_SQRTPI * t * exp(-t * t);
    }
    _Xv[N - 1] = 1.0;
}

////////////////////////////////////////////////////////////////////

double GaussianGeometry::density(double r) const
{
    double r2 = r * r;
    double sigma2 = _sigma * _sigma;
    return _rho0 * exp(-0.5 * r2 / sigma2);
}

////////////////////////////////////////////////////////////////////

double GaussianGeometry::randomRadius() const
{
    return random()->cdfLinLin(_rv, _Xv);
}

////////////////////////////////////////////////////////////////////

double GaussianGeometry::Sigmar() const
{
    return 1.0 / (4.0 * M_PI * _sigma * _sigma);
}

////////////////////////////////////////////////////////////////////
