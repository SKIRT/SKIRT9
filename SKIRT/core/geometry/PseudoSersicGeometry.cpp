/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "PseudoSersicGeometry.hpp"
#include "NR.hpp"
#include "Random.hpp"
#include "SpecialFunctions.hpp"

//////////////////////////////////////////////////////////////////////

void PseudoSersicGeometry::setupSelfBefore()
{
    SpheGeometry::setupSelfBefore();

    // calculate cached values
    _bn = 2.0 * _n - 1.0 / 3.0 + 4.0 / 405.0 / _n + 46.0 / 25515.0 / (_n * _n) + 131.0 / 1148175.0 / (_n * _n * _n);
    _rhon =
        pow(_bn, 2.0 * _n + 0.5) / (4.0 * M_PI * _n * SpecialFunctions::gamma(2.0 * _n + 0.5) * _reff * _reff * _reff);

    // grid with values of the cumulative mass
    int N = 201;
    double logsmin = -5.0;
    double logsmax = 5.0;
    double dlogs = (logsmax - logsmin) / (N - 1.0);
    _rv.resize(N);
    _Xv.resize(N);
    for (int i = 0; i < N; i++)
    {
        double logs = logsmin + i * dlogs;
        double s = pow(10.0, logs);
        _rv[i] = s * _reff;
        _Xv[i] = SpecialFunctions::incompleteGamma(2.0 * _n + 0.5, _bn * pow(s, 1.0 / _n));
    }
    _Xv[0] = 0.0;
    _Xv[N - 1] = 1.0;
}

////////////////////////////////////////////////////////////////////

double PseudoSersicGeometry::density(double r) const
{
    double s = r / _reff;
    double alpha = (2.0 * _n - 1.0) / (2.0 * _n);
    return _rhon * exp(-_bn * pow(s, 1.0 / _n)) * pow(s, -alpha);
}

//////////////////////////////////////////////////////////////////////

double PseudoSersicGeometry::randomRadius() const
{
    return random()->cdfLinLin(_rv, _Xv);
}

//////////////////////////////////////////////////////////////////////

double PseudoSersicGeometry::Sigmar() const
{
    return _n * _rhon * sqrt(M_PI) * _reff / sqrt(_bn);
}

///////////////////////////////////////////////////////////////////////////////
