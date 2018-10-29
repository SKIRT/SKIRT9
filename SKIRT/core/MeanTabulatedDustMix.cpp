/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "MeanTabulatedDustMix.hpp"
#include "FatalError.hpp"
#include "NR.hpp"

//////////////////////////////////////////////////////////////////////

void MeanTabulatedDustMix::setupSelfBefore()
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

MaterialMix::ScatteringMode MeanTabulatedDustMix::scatteringMode() const
{
    return ScatteringMode::HenyeyGreenstein;
}

//////////////////////////////////////////////////////////////////////

double MeanTabulatedDustMix::mass() const
{
    return _mu;
}

//////////////////////////////////////////////////////////////////////

double MeanTabulatedDustMix::sectionAbsSelf(double lambda) const
{
    size_t i = std::lower_bound(begin(_lambdav), end(_lambdav), lambda) - begin(_lambdav);
    if (i==0) return _sectionAbsv[0];
    if (i==_num) return _sectionAbsv[_num-1];
    return NR::interpolateLogLog(lambda, _lambdav[i-1], _lambdav[i], _sectionAbsv[i-1], _sectionAbsv[i]);
}

//////////////////////////////////////////////////////////////////////

double MeanTabulatedDustMix::sectionScaSelf(double lambda) const
{
    size_t i = std::lower_bound(begin(_lambdav), end(_lambdav), lambda) - begin(_lambdav);
    if (i==0) return _sectionScav[0];
    if (i==_num) return _sectionScav[_num-1];
    return NR::interpolateLogLog(lambda, _lambdav[i-1], _lambdav[i], _sectionScav[i-1], _sectionScav[i]);
}

//////////////////////////////////////////////////////////////////////

double MeanTabulatedDustMix::asymmpar(double lambda) const
{
    size_t i = std::lower_bound(begin(_lambdav), end(_lambdav), lambda) - begin(_lambdav);
    if (i==0) return _asymmparv[0];
    if (i==_num) return _asymmparv[_num-1];
    return NR::interpolateLinLin(lambda, _lambdav[i-1], _lambdav[i], _asymmparv[i-1], _asymmparv[i]);
}

//////////////////////////////////////////////////////////////////////
