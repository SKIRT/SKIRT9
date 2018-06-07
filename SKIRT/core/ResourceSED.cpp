/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "ResourceSED.hpp"
#include "Random.hpp"
#include "WavelengthRangeInterface.hpp"

//////////////////////////////////////////////////////////////////////

void ResourceSED::setupSelfBefore()
{
    SED::setupSelfBefore();

    _table.open(this, resourceName(), "lambda(m)", "Llambda(W/m)");
    _Ltot = _table.cdf(_lambdav, _cdfv, 10000, interface<WavelengthRangeInterface>()->wavelengthRange());
}

//////////////////////////////////////////////////////////////////////

double ResourceSED::specificLuminosity(double wavelength) const
{
    return _table[wavelength] / _Ltot;
}

//////////////////////////////////////////////////////////////////////

double ResourceSED::integratedLuminosity(const Range& wavelengthRange) const
{
    Array lambdav, cdfv;  // the contents of these arrays is not used, so this could be optimized if needed
    return _table.cdf(lambdav, cdfv, 1, wavelengthRange) / _Ltot;
}

//////////////////////////////////////////////////////////////////////

double ResourceSED::generateWavelength() const
{
    return random()->cdf(_lambdav, _cdfv);
}

//////////////////////////////////////////////////////////////////////
