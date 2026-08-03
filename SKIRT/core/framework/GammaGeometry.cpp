/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "GammaGeometry.hpp"
#include "Random.hpp"

//////////////////////////////////////////////////////////////////////

void GammaGeometry::setupSelfBefore()
{
    SpheGeometry::setupSelfBefore();

    // calculate cached values
    _rho0 = (3.0 - _gamma) / (4.0 * M_PI) / pow(_b, 3);
}

////////////////////////////////////////////////////////////////////

double GammaGeometry::density(double r) const
{
    return _rho0 * pow(r / _b, -_gamma) * pow(1.0 + r / _b, _gamma - 4.0);
}

//////////////////////////////////////////////////////////////////////

double GammaGeometry::randomRadius() const
{
    double X = random()->uniform();
    double t = pow(X, 1.0 / (3.0 - _gamma));
    return _b * t / (1.0 - t);
}

//////////////////////////////////////////////////////////////////////

double GammaGeometry::Sigmar() const
{
    if (_gamma < 1.0)
        return 1.0 / (2.0 * M_PI * _b * _b * (1.0 - _gamma) * (2.0 - _gamma));
    else
        return std::numeric_limits<double>::infinity();
}

//////////////////////////////////////////////////////////////////////
