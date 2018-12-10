/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "TabulatedDustMix.hpp"
#include "FatalError.hpp"
#include "NR.hpp"

//////////////////////////////////////////////////////////////////////

void TabulatedDustMix::setupSelfBefore()
{
    DustMix::setupSelfBefore();

    // obtain the wavelengths and properties
    Array kappaextv, albedov;
    _mu = getDustProperties(_lambdav, kappaextv, albedov, _asymmparv);

    // verify the number of points
    _num = _lambdav.size();
    if (_num < 1) throw FATALERROR("Dust properties must be tabulated for at least one wavelength");

    // calculate the cross sections
    _sectionScav = _mu * kappaextv * albedov;
    _sectionAbsv = _mu * kappaextv * (1.-albedov);
}

////////////////////////////////////////////////////////////////////

MaterialMix::ScatteringMode TabulatedDustMix::scatteringMode() const
{
    return ScatteringMode::HenyeyGreenstein;
}

//////////////////////////////////////////////////////////////////////

double TabulatedDustMix::mass() const
{
    return _mu;
}

//////////////////////////////////////////////////////////////////////

double TabulatedDustMix::sectionAbsSelf(double lambda) const
{
    return NR::clampedValue<NR::interpolateLogLog>(lambda, _lambdav, _sectionAbsv);
}

//////////////////////////////////////////////////////////////////////

double TabulatedDustMix::sectionScaSelf(double lambda) const
{
    return NR::clampedValue<NR::interpolateLogLog>(lambda, _lambdav, _sectionScav);
}

//////////////////////////////////////////////////////////////////////

double TabulatedDustMix::asymmpar(double lambda) const
{
    return NR::clampedValue<NR::interpolateLogLin>(lambda, _lambdav, _asymmparv);
}

//////////////////////////////////////////////////////////////////////
