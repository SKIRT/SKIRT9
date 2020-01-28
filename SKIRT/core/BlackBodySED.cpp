/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "BlackBodySED.hpp"
#include "PlanckFunction.hpp"
#include "Random.hpp"

//////////////////////////////////////////////////////////////////////

void BlackBodySED::setupSelfBefore()
{
    SED::setupSelfBefore();

    _planck = new PlanckFunction(temperature());
    _Ltot = _planck->cdf(_lambdav, _pv, _Pv, normalizationWavelengthRange());
}

//////////////////////////////////////////////////////////////////////

BlackBodySED::~BlackBodySED()
{
    delete _planck;
}

//////////////////////////////////////////////////////////////////////

Range BlackBodySED::intrinsicWavelengthRange() const
{
    return Range(std::numeric_limits<double>::denorm_min(), std::numeric_limits<double>::max());
}

//////////////////////////////////////////////////////////////////////

double BlackBodySED::specificLuminosity(double wavelength) const
{
    return _planck->value(wavelength) / _Ltot;
}

//////////////////////////////////////////////////////////////////////

void BlackBodySED::specificLuminosityArray(Array& lambdav, Array& pv, const Range& wavelengthRange) const
{
    Array Pv;  // the contents of this array is not used, so this could be optimized if needed
    double Ltot = _planck->cdf(lambdav, pv, Pv, wavelengthRange);
    pv *= (Ltot / _Ltot);
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
