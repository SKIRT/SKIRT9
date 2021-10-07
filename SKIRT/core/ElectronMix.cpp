/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "ElectronMix.hpp"
#include "Constants.hpp"
#include "MaterialState.hpp"
#include "PhotonPacket.hpp"
#include "Random.hpp"

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

bool ElectronMix::hasPolarizedScattering() const
{
    return includePolarization();
}

////////////////////////////////////////////////////////////////////

vector<StateVariable> ElectronMix::specificStateVariableInfo() const
{
    return vector<StateVariable>{StateVariable::numberDensity()};
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

double ElectronMix::opacityAbs(double /*lambda*/, const MaterialState* /*state*/, const PhotonPacket* /*pp*/) const
{
    return 0.;
}

////////////////////////////////////////////////////////////////////

double ElectronMix::opacitySca(double /*lambda*/, const MaterialState* state, const PhotonPacket* /*pp*/) const
{
    return state->numberDensity() * Constants::sigmaThomson();
}

////////////////////////////////////////////////////////////////////

double ElectronMix::opacityExt(double /*lambda*/, const MaterialState* state, const PhotonPacket* /*pp*/) const
{
    return state->numberDensity() * Constants::sigmaThomson();
}

////////////////////////////////////////////////////////////////////

void ElectronMix::peeloffScattering(double& I, double& Q, double& U, double& V, double& /*lambda*/, double w,
                                    Direction bfkobs, Direction bfky, const MaterialState* /*state*/,
                                    const PhotonPacket* pp) const
{
    _dpf.peeloffScattering(I, Q, U, V, w, pp->direction(), bfkobs, bfky, pp);
}

////////////////////////////////////////////////////////////////////

void ElectronMix::performScattering(double lambda, const MaterialState* state, PhotonPacket* pp) const
{
    // determine the new propagation direction, and if required, update the polarization state of the photon packet
    Direction bfknew = _dpf.performScattering(pp->direction(), pp);

    // execute the scattering event in the photon packet
    pp->scatter(bfknew, state->bulkVelocity(), lambda);
}

////////////////////////////////////////////////////////////////////
