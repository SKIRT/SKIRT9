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
    // get scaled energy with different definition
    double e = 2. * x;

    // draw random number from the distribution of the inverse Compton factor
    double r;
    while (true)
    {
        double xi1 = _random->uniform();
        double xi2 = _random->uniform();
        double xi3 = _random->uniform();

        if (xi1 <= 27. / (2. * e + 29.))
        {
            r = (e + 1.) / (e * xi2 + 1.);
            if (xi3 <= (std::pow((e + 2. - 2. * r) / e, 2) + 1.) / 2.) break;
        }
        else
        {
            r = e * xi2 + 1;
            if (xi3 <= 6.75 * (r - 1) * (r - 1) / (r * r * r)) break;
        }
    }

    // convert from inverse Compton factor to cosine theta
    return 1 - (r - 1) / x;
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
