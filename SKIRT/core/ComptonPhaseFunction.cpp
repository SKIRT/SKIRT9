/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "ComptonPhaseFunction.hpp"
#include "Constants.hpp"
#include "NR.hpp"
#include "Random.hpp"

////////////////////////////////////////////////////////////////////

namespace
{
    // discretization of the phase function over scattering angle: theta from 0 to pi, index t
    constexpr int numTheta = 361;
    constexpr int maxTheta = numTheta - 1;
    constexpr double deltaTheta = M_PI / maxTheta;

    // returns the photon energy scaled to the electron rest energy: h nu / m_e c^2
    double scaledEnergy(double lambda) { return (Constants::h() / Constants::Melectron() / Constants::c()) / lambda; }

    // returns the Compton scattering cross section for a given scaled energy
    double comptonSection(double x)
    {
        double x2 = x * x;
        double x3 = x2 * x;
        double x1 = 1. + x;
        double x12 = 1. + 2. * x;
        return 0.75 * (x1 / (x12 * x12) + 2. / x2 + (1. / (2. * x) - x1 / x3) * log(x12));
    }

    // returns the inverse Compton factor for a given scaled energy and scattering angle cosine
    double inverseComptonfactor(double x, double costheta) { return 1 + x * (1 - costheta); }

    // returns the Compton factor for a given scaled energy and scattering angle cosine
    double comptonFactor(double x, double costheta) { return 1. / inverseComptonfactor(x, costheta); }
}

////////////////////////////////////////////////////////////////////

void ComptonPhaseFunction::initialize(Random* random)
{
    // cache random number generator
    _random = random;

    // construct a theta grid and precalculate values used in generateCosineFromPhaseFunction()
    // to accelerate construction of the cumulative phase function distribution
    _costhetav.resize(numTheta);
    _sinthetav.resize(numTheta);
    _sin2thetav.resize(numTheta);
    for (int t = 0; t != numTheta; ++t)
    {
        double theta = t * deltaTheta;
        _costhetav[t] = cos(theta);
        _sinthetav[t] = sin(theta);
        _sin2thetav[t] = _sinthetav[t] * _sinthetav[t];
    }
}

////////////////////////////////////////////////////////////////////

double ComptonPhaseFunction::sectionSca(double lambda) const
{
    double x = scaledEnergy(lambda);
    return comptonSection(x);
}

////////////////////////////////////////////////////////////////////

double ComptonPhaseFunction::phaseFunctionValueForCosine(double x, double costheta) const
{
    double C = comptonFactor(x, costheta);
    double sin2theta = (1 - costheta) * (1 + costheta);
    double phase = C * C * C + C - C * C * sin2theta;
    return 0.75 / comptonSection(x) * phase;
}

////////////////////////////////////////////////////////////////////

double ComptonPhaseFunction::generateCosineFromPhaseFunction(double x) const
{
    // construct the normalized cumulative phase function distribution for this x
    Array thetaXv;
    NR::cdf(thetaXv, maxTheta, [this, x](int t) {
        t += 1;
        double C = comptonFactor(x, _costhetav[t]);
        double phase = C * C * C + C - C * C * _sin2thetav[t];
        return phase * _sinthetav[t];
    });

    // draw a random cosine from this distribution
    return _random->cdfLinLin(_costhetav, thetaXv);
}

////////////////////////////////////////////////////////////////////

void ComptonPhaseFunction::peeloffScattering(double& I, double& lambda, Direction bfk, Direction bfkobs) const
{
    double x = scaledEnergy(lambda);

    // calculate the value of the phase function
    double costheta = Vec::dot(bfk, bfkobs);
    double value = phaseFunctionValueForCosine(x, costheta);

    // accumulate the weighted sum in the intensity
    I += value;

    // adjust the wavelength
    lambda *= inverseComptonfactor(x, costheta);
}

////////////////////////////////////////////////////////////////////

Direction ComptonPhaseFunction::performScattering(double& lambda, Direction bfk) const
{
    double x = scaledEnergy(lambda);

    // sample a scattering angle from the phase function
    double costheta = generateCosineFromPhaseFunction(x);

    // adjust the wavelength
    lambda *= inverseComptonfactor(x, costheta);

    // determine the new propagation direction
    return _random->direction(bfk, costheta);
}

////////////////////////////////////////////////////////////////////
