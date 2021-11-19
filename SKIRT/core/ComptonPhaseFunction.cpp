/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "ComptonPhaseFunction.hpp"
#include "Random.hpp"

////////////////////////////////////////////////////////////////////

void ComptonPhaseFunction::initialize(Random* random)
{
    // cache random number generator
    _random = random;
}

////////////////////////////////////////////////////////////////////

double ComptonPhaseFunction::sectionSca(double /*lambda*/) const
{
    return 1.;
}

////////////////////////////////////////////////////////////////////

double ComptonPhaseFunction::phaseFunctionValueForCosine(double costheta) const
{
    return 0.75 * (costheta * costheta + 1.);
}

////////////////////////////////////////////////////////////////////

double ComptonPhaseFunction::generateCosineFromPhaseFunction() const
{
    double X = _random->uniform();
    double p = cbrt(4. * X - 2. + sqrt(16. * X * (X - 1.) + 5.));
    return p - 1. / p;
}

////////////////////////////////////////////////////////////////////

void ComptonPhaseFunction::peeloffScattering(double& I, double& /*lambda*/, double w, Direction bfk,
                                             Direction bfkobs) const
{
    // calculate the value of the material-specific phase function
    double costheta = Vec::dot(bfk, bfkobs);
    double value = phaseFunctionValueForCosine(costheta);

    // accumulate the weighted sum in the intensity
    I += w * value;
}

////////////////////////////////////////////////////////////////////

Direction ComptonPhaseFunction::performScattering(double& /*lambda*/, Direction bfk) const
{
    // sample a scattering angle from the dipole phase function
    double costheta = generateCosineFromPhaseFunction();

    // determine the new propagation direction
    return _random->direction(bfk, costheta);
}

////////////////////////////////////////////////////////////////////
