/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "EinastoGeometry.hpp"
#include "Random.hpp"
#include "SpecialFunctions.hpp"

//////////////////////////////////////////////////////////////////////

void EinastoGeometry::setupSelfBefore()
{
    SpheGeometry::setupSelfBefore();

    // calculate cached values
    _dn = 3.0 * _n - 1.0 / 3.0 + 8.0 / 1215.0 / _n + 184.0 / 229635.0 / _n / _n + 1048.0 / 31000725.0 / pow(_n, 3);
    _rhos = 1.0 / (4.0 * M_PI) / pow(_rs, 3) / _n / SpecialFunctions::gamma(3.0 * _n) / exp(_dn) * pow(_dn, 3.0 * _n);

    // grid with values of the cumulative mass
    int N = 401;
    double logsmin = -4.0;
    double logsmax = 4.0;
    double dlogs = (logsmax - logsmin) / (N - 1.0);
    _rv.resize(N);
    _Xv.resize(N);
    for (int i = 0; i < N; i++)
    {
        double logs = logsmin + i * dlogs;
        double s = pow(10.0, logs);
        _rv[i] = s * _rs;
        _Xv[i] = SpecialFunctions::incompleteGamma(3.0 * _n, _dn * pow(s, 1.0 / _n));
    }
    _Xv[0] = 0.0;
    _Xv[N - 1] = 1.0;
}

//////////////////////////////////////////////////////////////////////

double EinastoGeometry::density(double r) const
{
    return _rhos * exp(-_dn * (pow(r / _rs, 1.0 / _n) - 1.0));
}

//////////////////////////////////////////////////////////////////////

double EinastoGeometry::randomRadius() const
{
    return random()->cdfLinLin(_rv, _Xv);
}

//////////////////////////////////////////////////////////////////////

double EinastoGeometry::Sigmar() const
{
    return _n * SpecialFunctions::gamma(_n) * _rhos * exp(_dn) * _rs / pow(_dn, _n);
}

///////////////////////////////////////////////////////////////////////////////
