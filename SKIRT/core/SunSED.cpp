/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "SunSED.hpp"
#include "Random.hpp"

//////////////////////////////////////////////////////////////////////

void SunSED::setupSelfBefore()
{
    SED::setupSelfBefore();

    _table.open(this, "SunSED", "lambda(m)", "Llambda(W/m)");
    _Ltot = _table.cdf(_lambdav, _cdfv, 200, minWavelength(), maxWavelength());
}

//////////////////////////////////////////////////////////////////////

double SunSED::specificLuminosity(double wavelength) const
{
    return _table[wavelength] / _Ltot;
}

//////////////////////////////////////////////////////////////////////

double SunSED::integratedLuminosity(double minWavelength, double maxWavelength) const
{
    Array lambdav, cdfv;  // the contents of these arrays is not used, so this could be optimized if needed
    return _table.cdf(lambdav, cdfv, 1, minWavelength, maxWavelength) / _Ltot;
}

//////////////////////////////////////////////////////////////////////

double SunSED::generateWavelength() const
{
    return random()->cdf(_lambdav, _cdfv);
}

//////////////////////////////////////////////////////////////////////
