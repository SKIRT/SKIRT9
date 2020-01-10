/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "SersicGeometry.hpp"
#include "FatalError.hpp"
#include "Random.hpp"
#include "SersicFunction.hpp"
#include "SpecialFunctions.hpp"

//////////////////////////////////////////////////////////////////////

SersicGeometry::~SersicGeometry()
{
    delete _sersicfunction;
}

//////////////////////////////////////////////////////////////////////

void SersicGeometry::setupSelfBefore()
{
    SpheGeometry::setupSelfBefore();

    // calculate cached values
    _rho0 = 1.0 / (_reff * _reff * _reff);
    _b = 2.0 * _n - 1.0 / 3.0 + 4.0 / 405.0 / _n + 46.0 / 25515.0 / (_n * _n) + 131.0 / 1148175.0 / (_n * _n * _n);
    _sersicfunction = new SersicFunction(_n);
}

//////////////////////////////////////////////////////////////////////

double SersicGeometry::density(double r) const
{
    double s = r / _reff;
    return _rho0 * (*_sersicfunction)(s);
}

//////////////////////////////////////////////////////////////////////

double SersicGeometry::randomRadius() const
{
    double X = random()->uniform();
    return _reff * _sersicfunction->inverseMass(X);
}

//////////////////////////////////////////////////////////////////////

double SersicGeometry::Sigmar() const
{
    return 1.0 / (_reff * _reff) * pow(_b, 2.0 * _n) / (2.0 * M_PI * SpecialFunctions::gamma(2.0 * _n + 1.0));
}

//////////////////////////////////////////////////////////////////////
