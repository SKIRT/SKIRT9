/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "FamilySED.hpp"
#include "Random.hpp"
#include "SEDFamily.hpp"

//////////////////////////////////////////////////////////////////////

void FamilySED::setupSelfBefore()
{
    SED::setupSelfBefore();

    _family = getFamilyAndParameters(_parameters);
    _Ltot = _family->cdf(_lambdav, _pv, _Pv, normalizationWavelengthRange(), _parameters);
}

//////////////////////////////////////////////////////////////////////

Range FamilySED::intrinsicWavelengthRange() const
{
    return _family->intrinsicWavelengthRange();
}

//////////////////////////////////////////////////////////////////////

double FamilySED::specificLuminosity(double wavelength) const
{
    return _family->specificLuminosity(wavelength, _parameters) / _Ltot;
}

//////////////////////////////////////////////////////////////////////

void FamilySED::specificLuminosityArray(Array& lambdav, Array& pv, const Range& wavelengthRange) const
{
    Array Pv;  // the contents of this array is not used, so this could be optimized if needed
    double Ltot = _family->cdf(lambdav, pv, Pv, wavelengthRange, _parameters);
    pv *= (Ltot / _Ltot);
}

//////////////////////////////////////////////////////////////////////

double FamilySED::integratedLuminosity(const Range& wavelengthRange) const
{
    Array lambdav, pv, Pv;  // the contents of these arrays is not used, so this could be optimized if needed
    return _family->cdf(lambdav, pv, Pv, wavelengthRange, _parameters) / _Ltot;
}

//////////////////////////////////////////////////////////////////////

double FamilySED::generateWavelength() const
{
    return random()->cdfLogLog(_lambdav, _pv, _Pv);
}

//////////////////////////////////////////////////////////////////////
