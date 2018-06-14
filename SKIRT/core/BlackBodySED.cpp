/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "BlackBodySED.hpp"
#include "Random.hpp"
#include "WavelengthRangeInterface.hpp"
#include "PlanckFunction.hpp"

//////////////////////////////////////////////////////////////////////

void BlackBodySED::setupSelfBefore()
{
    SED::setupSelfBefore();

    _planck = new PlanckFunction(temperature());
    _Ltot = _planck->cdf(_lambdav, _pv, _Pv, interface<WavelengthRangeInterface>()->wavelengthRange());
}

//////////////////////////////////////////////////////////////////////

BlackBodySED::~BlackBodySED()
{
    delete _planck;
}

//////////////////////////////////////////////////////////////////////

double BlackBodySED::specificLuminosity(double wavelength) const
{
    return _planck->value(wavelength) / _Ltot;
}

//////////////////////////////////////////////////////////////////////

double BlackBodySED::integratedLuminosity(const Range& wavelengthRange) const
{
    Array lambdav, pv, Pv;  // the contents of these arrays is not used, so this could be optimized if needed
    return _planck->cdf(lambdav, pv, Pv, wavelengthRange) / _Ltot;
}

//////////////////////////////////////////////////////////////////////

double BlackBodySED::generateWavelength() const
{
    return random()->cdfLogLog(_lambdav, _pv, _Pv);
}

//////////////////////////////////////////////////////////////////////
