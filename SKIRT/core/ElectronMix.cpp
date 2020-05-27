/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "ElectronMix.hpp"
#include "Constants.hpp"
#include "Random.hpp"
#include "StokesVector.hpp"

////////////////////////////////////////////////////////////////////

void ElectronMix::setupSelfBefore()
{
    MaterialMix::setupSelfBefore();

    _dpf.initialize(random(), includePolarization());
}

////////////////////////////////////////////////////////////////////

MaterialMix::MaterialType ElectronMix::materialType() const
{
    return MaterialType::Electrons;
}

////////////////////////////////////////////////////////////////////

MaterialMix::ScatteringMode ElectronMix::scatteringMode() const
{
    return includePolarization() ? ScatteringMode::SphericalPolarization : ScatteringMode::MaterialPhaseFunction;
}

////////////////////////////////////////////////////////////////////

double ElectronMix::mass() const
{
    return Constants::Melectron();
}

////////////////////////////////////////////////////////////////////

double ElectronMix::sectionAbs(double /*lambda*/) const
{
    return 0.;
}

////////////////////////////////////////////////////////////////////

double ElectronMix::sectionSca(double /*lambda*/) const
{
    return Constants::sigmaThomson();
}

////////////////////////////////////////////////////////////////////

double ElectronMix::sectionExt(double /*lambda*/) const
{
    return Constants::sigmaThomson();
}

////////////////////////////////////////////////////////////////////

double ElectronMix::phaseFunctionValueForCosine(double /*lambda*/, double costheta) const
{
    return _dpf.phaseFunctionValueForCosine(costheta);
}

////////////////////////////////////////////////////////////////////

double ElectronMix::generateCosineFromPhaseFunction(double /*lambda*/) const
{
    return _dpf.generateCosineFromPhaseFunction();
}

////////////////////////////////////////////////////////////////////

double ElectronMix::phaseFunctionValue(double /*lambda*/, double theta, double phi, const StokesVector* sv) const
{
    return _dpf.phaseFunctionValue(theta, phi, sv);
}

////////////////////////////////////////////////////////////////////

std::pair<double, double> ElectronMix::generateAnglesFromPhaseFunction(double /*lambda*/, const StokesVector* sv) const
{
    return _dpf.generateAnglesFromPhaseFunction(sv);
}

////////////////////////////////////////////////////////////////////

void ElectronMix::applyMueller(double /*lambda*/, double theta, StokesVector* sv) const
{
    _dpf.applyMueller(theta, sv);
}

////////////////////////////////////////////////////////////////////

double ElectronMix::equilibriumTemperature(const Array& /*Jv*/) const
{
    return 0.;
}

////////////////////////////////////////////////////////////////////

Array ElectronMix::emissivity(const Array& /*Jv*/) const
{
    return Array();
}

////////////////////////////////////////////////////////////////////
