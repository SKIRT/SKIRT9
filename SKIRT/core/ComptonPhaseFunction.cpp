/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "ComptonPhaseFunction.hpp"
#include "Constants.hpp"
#include "Random.hpp"

////////////////////////////////////////////////////////////////////

namespace
{
    // returns the photon energy scaled to the electron rest energy: h nu / m_e c^2
    double scaledEnergy(double lambda) { return (Constants::h() / Constants::Melectron() / Constants::c()) / lambda; }

    // returns the inverse Compton factor for a given scaled energy and scattering angle cosine
    double inverseComptonfactor(double x, double costheta) { return 1 + x * (1 - costheta); }

    // returns the Compton factor for a given scaled energy and scattering angle cosine
    ///double straightComptonfactor(double x, double costheta) { return 1. / inverseComptonfactor(x, costheta); }
}

////////////////////////////////////////////////////////////////////

void ComptonPhaseFunction::initialize(Random* random)
{
    // cache random number generator
    _random = random;
}

////////////////////////////////////////////////////////////////////

double ComptonPhaseFunction::sectionSca(double lambda) const
{
    double x = scaledEnergy(lambda);
    double x2 = x * x;
    double x3 = x2 * x;
    double x1 = 1. + x;
    double x12 = 1. + 2. * x;
    return 0.75 * (x1 / (x12 * x12) + 2. / x2 + (1. / (2. * x) - x1 / x3) * log(x12));
}

////////////////////////////////////////////////////////////////////

double ComptonPhaseFunction::phaseFunctionValueForCosine(double /*x*/, double costheta) const
{
    return 0.75 * (costheta * costheta + 1.);
}

////////////////////////////////////////////////////////////////////

double ComptonPhaseFunction::generateCosineFromPhaseFunction(double /*x*/) const
{
    double X = _random->uniform();
    double p = cbrt(4. * X - 2. + sqrt(16. * X * (X - 1.) + 5.));
    return p - 1. / p;
}

////////////////////////////////////////////////////////////////////

void ComptonPhaseFunction::peeloffScattering(double& I, double& lambda, double w, Direction bfk, Direction bfkobs) const
{
    double x = scaledEnergy(lambda);

    // calculate the value of the phase function
    double costheta = Vec::dot(bfk, bfkobs);
    double value = phaseFunctionValueForCosine(x, costheta);

    // accumulate the weighted sum in the intensity
    I += w * value;

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
