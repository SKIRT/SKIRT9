/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "DustMix.hpp"
#include "Configuration.hpp"
#include "Log.hpp"
#include "NR.hpp"
#include "Random.hpp"
#include "StokesVector.hpp"
#include "StringUtils.hpp"

////////////////////////////////////////////////////////////////////

// built-in constants determining the resolution for discretizing the scattering angles
namespace
{
    // theta from 0 to pi, index t
    constexpr int numTheta = 361;
    constexpr int maxTheta = numTheta-1;
    constexpr double deltaTheta = M_PI/maxTheta;

    // phi from 0 to 2 pi, index f
    constexpr int numPhi = 361;
    constexpr int maxPhi = numPhi-1;
    constexpr double deltaPhi = 2*M_PI/(maxPhi-1);
}

////////////////////////////////////////////////////////////////////

void DustMix::setupSelfAfter()
{
    MaterialMix::setupSelfAfter();

    // get the scattering mode advertised by this dust mix
    auto mode = scatteringMode();

    // determine the parameters for a fine grid covering the wavelength range of the simulation in log space;
    // use integer multiples as logarithmic grid points so that the grid is stable for changing wavelength ranges
    const int numWavelengthsPerDex = 1024;
    Range range = find<Configuration>()->simulationWavelengthRange();
    int minLambdaSerial = std::floor(numWavelengthsPerDex*log10(range.min()));
    int maxLambdaSerial = std::ceil(numWavelengthsPerDex*log10(range.max()));
    _logLambdaFactor = numWavelengthsPerDex;
    _logLambdaOffset = minLambdaSerial;
    _maxLambdaIndex = maxLambdaSerial - minLambdaSerial;
    int numLambda = _maxLambdaIndex+1;

    // build a temporary wavelength grid corresponding to this scheme
    Array lambdav(numLambda);
    for (int ell=0; ell!=numLambda; ++ell) lambdav[ell] = pow(10., (ell+_logLambdaOffset+0.5)/_logLambdaFactor);

    // if needed, build a scattering angle grid
    if (mode == ScatteringMode::MaterialPhaseFunction || mode == ScatteringMode::SphericalPolarization)
    {
        _thetav.resize(numTheta);
        for (int t=0; t!=numTheta; ++t) _thetav[t] = t * deltaTheta;
    }

    // resize the optical property arrays and tables as needed
    _sigmaabsv.resize(numLambda);
    _sigmascav.resize(numLambda);
    _sigmaextv.resize(numLambda);
    _albedov.resize(numLambda);
    _asymmparv.resize(numLambda);
    if (mode == ScatteringMode::MaterialPhaseFunction || mode == ScatteringMode::SphericalPolarization)
    {
        _S11vv.resize(numLambda, numTheta);
        if (mode == ScatteringMode::SphericalPolarization)
        {
            _S12vv.resize(numLambda, numTheta);
            _S33vv.resize(numLambda, numTheta);
            _S34vv.resize(numLambda, numTheta);
        }
    }

    // obtain the optical properties from the subclass
    _mu = getOpticalProperties(lambdav, _thetav, _sigmaabsv, _sigmascav, _asymmparv, _S11vv, _S12vv, _S33vv, _S34vv);

    // calculate some derived basic optical properties
    for (int ell=0; ell!=numLambda; ++ell)
    {
        _sigmaextv[ell] = _sigmaabsv[ell] + _sigmascav[ell];
        _albedov[ell] = _sigmaextv[ell] > 0. ? _sigmascav[ell]/_sigmaextv[ell] : 0.;
    }

    // precalculate discretizations related to the scattering angles as needed
    if (mode == ScatteringMode::MaterialPhaseFunction || mode == ScatteringMode::SphericalPolarization)
    {
        // create a table with the normalized cumulative distribution of theta for each wavelength
        _thetaXvv.resize(numLambda,0);
        for (int ell=0; ell!=numLambda; ++ell)
        {
            NR::cdf(_thetaXvv[ell], maxTheta, [this,ell](int t){ return _S11vv(ell,t+1)*sin(_thetav[t+1]); });
        }

        // create a table with the phase function normalization factor for each wavelength
        _pfnormv.resize(numLambda);
        for (int ell=0; ell!=numLambda; ++ell)
        {
            double sum = 0.;
            for (int t=0; t!=numTheta; ++t)
            {
                sum += _S11vv(ell,t)*sin(_thetav[t])*deltaTheta;
            }
            _pfnormv[ell] = 2.0/sum;
        }

        // create tables listing phi, phi/(2 pi), sin(2 phi) and 1-cos(2 phi) for each phi index
        if (mode == ScatteringMode::SphericalPolarization)
        {
            _phiv.resize(numPhi);
            _phi1v.resize(numPhi);
            _phisv.resize(numPhi);
            _phicv.resize(numPhi);
            for (int f=0; f!=numPhi; ++f)
            {
                double phi = f * deltaPhi;
                _phiv[f] = phi;
                _phi1v[f] = phi/(2*M_PI);
                _phisv[f] = sin(2*phi);
                _phicv[f] = 1-cos(2*phi);
            }
        }
    }

    // precalculate information to accelerate solving the energy balance equation for the temperature;
    // this is relevant only if the simulation tracks the radiation field
    if (find<Configuration>()->hasPanRadiationField())
    {
        _calc.precalculate(this, lambdav, _sigmaabsv);
    }

    // give the subclass a chance to obtain additional precalculated information
    size_t allocatedBytes = initializeExtraProperties(lambdav);

    // calculate and log allocated memory size
    size_t allocatedSize = 0;
    allocatedSize += _thetav.size();
    allocatedSize += _sigmaabsv.size();
    allocatedSize += _sigmascav.size();
    allocatedSize += _sigmaextv.size();
    allocatedSize += _albedov.size();
    allocatedSize += _asymmparv.size();
    allocatedSize += _S11vv.size();
    allocatedSize += _S12vv.size();
    allocatedSize += _S33vv.size();
    allocatedSize += _S34vv.size();
    allocatedSize += _thetaXvv.size();
    allocatedSize += _pfnormv.size();
    allocatedSize += _phiv.size();
    allocatedSize += _phi1v.size();
    allocatedSize += _phisv.size();
    allocatedSize += _phicv.size();

    allocatedBytes += allocatedSize*sizeof(double) + _calc.allocatedBytes();
    find<Log>()->info(type() + " allocated " + StringUtils::toMemSizeString(allocatedBytes) + " of memory");
}

////////////////////////////////////////////////////////////////////

size_t DustMix::initializeExtraProperties(const Array& /*lambdav*/)
{
    return 0;
}

////////////////////////////////////////////////////////////////////

int DustMix::indexForLambda(double lambda) const
{
    int ell = static_cast<int>(log10(lambda) * _logLambdaFactor) - _logLambdaOffset;
    return max(0, min(ell, _maxLambdaIndex));
}

////////////////////////////////////////////////////////////////////

int DustMix::indexForTheta(double theta) const
{
    int t = 0.5 + theta / deltaTheta;
    return max(0, min(t, maxTheta));
}

////////////////////////////////////////////////////////////////////

MaterialMix::MaterialType DustMix::materialType() const
{
    return MaterialType::Dust;
}

////////////////////////////////////////////////////////////////////

double DustMix::mass() const
{
    return _mu;
}

////////////////////////////////////////////////////////////////////

double DustMix::sectionAbs(double lambda) const
{
    return _sigmaabsv[indexForLambda(lambda)];
}

////////////////////////////////////////////////////////////////////

double DustMix::sectionSca(double lambda) const
{
    return _sigmascav[indexForLambda(lambda)];
}

////////////////////////////////////////////////////////////////////

double DustMix::sectionExt(double lambda) const
{
    return _sigmaextv[indexForLambda(lambda)];
}

////////////////////////////////////////////////////////////////////

double DustMix::albedo(double lambda) const
{
    return _albedov[indexForLambda(lambda)];
}

////////////////////////////////////////////////////////////////////

double DustMix::asymmpar(double lambda) const
{
    return _asymmparv[indexForLambda(lambda)];
}

////////////////////////////////////////////////////////////////////

double DustMix::phaseFunctionValueForCosine(double lambda, double costheta) const
{
    int ell = indexForLambda(lambda);
    int t = indexForTheta(acos(costheta));
    return _pfnormv[ell] * _S11vv(ell,t);
}

////////////////////////////////////////////////////////////////////

double DustMix::generateCosineFromPhaseFunction(double lambda) const
{
    return cos(random()->cdfLinLin(_thetav, _thetaXvv[indexForLambda(lambda)]));
}

////////////////////////////////////////////////////////////////////

double DustMix::phaseFunctionValue(double lambda, double theta, double phi, const StokesVector* sv) const
{
    int ell = indexForLambda(lambda);
    int t = indexForTheta(theta);
    double polDegree = sv->linearPolarizationDegree();
    double polAngle = sv->polarizationAngle();
    return _pfnormv[ell] * (_S11vv(ell,t) + polDegree*_S12vv(ell,t)*cos(2.*(phi-polAngle)));
}

////////////////////////////////////////////////////////////////////

std::pair<double,double> DustMix::generateAnglesFromPhaseFunction(double lambda, const StokesVector* sv) const
{
    int ell = indexForLambda(lambda);

    // sample from the normalized cumulative distribution of theta for this wavelength
    double theta = random()->cdfLinLin(_thetav, _thetaXvv[ell]);
    int t = indexForTheta(theta);

    // construct and sample from the normalized cumulative distribution of phi for this wavelength and theta angle
    double polDegree = sv->linearPolarizationDegree();
    double polAngle = sv->polarizationAngle();
    double PF = polDegree * _S12vv(ell,t)/_S11vv(ell,t) / (4*M_PI);
    double cos2polAngle = cos(2*polAngle) * PF;
    double sin2polAngle = sin(2*polAngle) * PF;
    double phi = random()->cdfLinLin(_phiv, _phi1v + cos2polAngle*_phisv + sin2polAngle*_phicv);

    // return the result
    return std::make_pair(theta, phi);
}

////////////////////////////////////////////////////////////////////

void DustMix::applyMueller(double lambda, double theta, StokesVector* sv) const
{
    int ell = indexForLambda(lambda);
    int t = indexForTheta(theta);
    sv->applyMueller(_S11vv(ell,t), _S12vv(ell,t), _S33vv(ell,t), _S34vv(ell,t));
}

////////////////////////////////////////////////////////////////////

double DustMix::equilibriumTemperature(const Array& Jv) const
{
    return _calc.equilibriumTemperature(0, Jv);
}

////////////////////////////////////////////////////////////////////

Array DustMix::emissivity(const Array& Jv) const
{
    return _calc.emissivity(Jv);
}

////////////////////////////////////////////////////////////////////
