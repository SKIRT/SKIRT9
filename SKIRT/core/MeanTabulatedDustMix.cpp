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
    MaterialMix::setupSelfBefore();

    // obtain the wavelengths and properties
    Array kappaextv;
    _mu = getDustProperties(_lambdav, kappaextv, _albedov, _asymmparv);

    // verify the number of points
    _num = _lambdav.size();
    if (_num < 1) throw FATALERROR("Dust properties must be tabulated for at least one wavelength");

    // calculate the cross sections
    _sectionExtv = kappaextv * _mu;
    _sectionScav = _sectionExtv * _albedov;
    _sectionAbsv = _sectionExtv * (1.-_albedov);
}

////////////////////////////////////////////////////////////////////

MaterialMix::MaterialType MeanTabulatedDustMix::materialType() const
{
    return MaterialType::Dust;
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

double MeanTabulatedDustMix::sectionAbs(double lambda) const
{
    size_t i = std::lower_bound(begin(_lambdav), end(_lambdav), lambda) - begin(_lambdav);
    if (i==0) return _sectionAbsv[0];
    if (i==_num) return _sectionAbsv[_num-1];
    return NR::interpolateLogLog(lambda, _lambdav[i-1], _lambdav[i], _sectionAbsv[i-1], _sectionAbsv[i]);
}

//////////////////////////////////////////////////////////////////////

double MeanTabulatedDustMix::sectionSca(double lambda) const
{
    size_t i = std::lower_bound(begin(_lambdav), end(_lambdav), lambda) - begin(_lambdav);
    if (i==0) return _sectionScav[0];
    if (i==_num) return _sectionScav[_num-1];
    return NR::interpolateLogLog(lambda, _lambdav[i-1], _lambdav[i], _sectionScav[i-1], _sectionScav[i]);
}

//////////////////////////////////////////////////////////////////////

double MeanTabulatedDustMix::sectionExt(double lambda) const
{
    size_t i = std::lower_bound(begin(_lambdav), end(_lambdav), lambda) - begin(_lambdav);
    if (i==0) return _sectionExtv[0];
    if (i==_num) return _sectionExtv[_num-1];
    return NR::interpolateLogLog(lambda, _lambdav[i-1], _lambdav[i], _sectionExtv[i-1], _sectionExtv[i]);
}

//////////////////////////////////////////////////////////////////////

double MeanTabulatedDustMix::albedo(double lambda) const
{
    size_t i = std::lower_bound(begin(_lambdav), end(_lambdav), lambda) - begin(_lambdav);
    if (i==0) return _albedov[0];
    if (i==_num) return _albedov[_num-1];
    return NR::interpolateLogLog(lambda, _lambdav[i-1], _lambdav[i], _albedov[i-1], _albedov[i]);
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
