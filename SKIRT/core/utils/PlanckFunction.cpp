/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "PlanckFunction.hpp"
#include "Constants.hpp"
#include "NR.hpp"

////////////////////////////////////////////////////////////////////

PlanckFunction::PlanckFunction(double T) : _T(T)
{
    // precompute some front factors
    double h = Constants::h();
    double c = Constants::c();
    double k = Constants::k();
    _f1 = h * c / (k * _T);
    _f2 = 2.0 * h * c * c;
}

////////////////////////////////////////////////////////////////////

double PlanckFunction::value(double lambda) const
{
    return _f2 / pow(lambda, 5) / (exp(_f1 / lambda) - 1.0);
}

////////////////////////////////////////////////////////////////////

double PlanckFunction::cdf(Array& lambdav, Array& pv, Array& Pv, Range lambdaRange) const
{
    // build an appropriate grid
    size_t n = max(static_cast<size_t>(100), static_cast<size_t>(1000. * log10(lambdaRange.max() / lambdaRange.min())));
    NR::buildLogGrid(lambdav, lambdaRange.min(), lambdaRange.max(), n);

    // calculate the tabulated probability densities
    pv.resize(n + 1);
    for (size_t i = 0; i <= n; ++i) pv[i] = value(lambdav[i]);

    // perform the rest of the operation in an implementation function using log-log interpolation
    return NR::cdf2(true, lambdav, pv, Pv);
}

////////////////////////////////////////////////////////////////////
