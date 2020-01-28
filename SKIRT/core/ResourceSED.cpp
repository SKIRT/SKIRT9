/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "ResourceSED.hpp"
#include "Random.hpp"

//////////////////////////////////////////////////////////////////////

void ResourceSED::setupSelfBefore()
{
    SED::setupSelfBefore();

    _table.open(this, resourceName(), "lambda(m)", "Llambda(W/m)", false);
    _Ltot = _table.cdf(_lambdav, _pv, _Pv, normalizationWavelengthRange());
}

//////////////////////////////////////////////////////////////////////

Range ResourceSED::intrinsicWavelengthRange() const
{
    return _table.axisRange<0>();
}

//////////////////////////////////////////////////////////////////////

double ResourceSED::specificLuminosity(double wavelength) const
{
    return _table[wavelength] / _Ltot;
}

//////////////////////////////////////////////////////////////////////

void ResourceSED::specificLuminosityArray(Array& lambdav, Array& pv, const Range& wavelengthRange) const
{
    Array Pv;  // the contents of this array is not used, so this could be optimized if needed
    double Ltot = _table.cdf(lambdav, pv, Pv, wavelengthRange);
    pv *= (Ltot / _Ltot);
}

//////////////////////////////////////////////////////////////////////

double ResourceSED::integratedLuminosity(const Range& wavelengthRange) const
{
    Array lambdav, pv, Pv;  // the contents of these arrays is not used, so this could be optimized if needed
    return _table.cdf(lambdav, pv, Pv, wavelengthRange) / _Ltot;
}

//////////////////////////////////////////////////////////////////////

double ResourceSED::generateWavelength() const
{
    return random()->cdfLogLog(_lambdav, _pv, _Pv);
}

//////////////////////////////////////////////////////////////////////
