/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "ElectronMix.hpp"
#include "Configuration.hpp"
#include "Constants.hpp"
#include "Log.hpp"
#include "MaterialState.hpp"
#include "PhotonPacket.hpp"
#include "Random.hpp"

////////////////////////////////////////////////////////////////////

namespace
{
    // transition wavelength from Compton to Thomson scattering
    const double comptonWL = 10e-9;  // 10 nm
}

////////////////////////////////////////////////////////////////////

void ElectronMix::setupSelfBefore()
{
    MaterialMix::setupSelfBefore();

    // determine whether we need Compton scattering because the simulation may have short wavelengths
    auto range = find<Configuration>()->simulationWavelengthRange();
    _hasCompton = range.min() < comptonWL;

    // disable Compton scattering if polarization is requested
    if (_hasCompton && includePolarization())
    {
        _hasCompton = false;
        find<Log>()->warning("Compton scattering disabled because the implementation does not support polarization");
    }

    // initialize our phase function helpers
    _dpf.initialize(random(), includePolarization());
    if (_hasCompton) _cpf.initialize(random());
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

double ElectronMix::sectionSca(double lambda) const
{
    double sigma = Constants::sigmaThomson();
    if (_hasCompton && lambda < comptonWL) sigma *= _cpf.sectionSca(lambda);
    return sigma;
}

////////////////////////////////////////////////////////////////////

double ElectronMix::sectionExt(double lambda) const
{
    return sectionSca(lambda);
}

////////////////////////////////////////////////////////////////////

double ElectronMix::opacityAbs(double /*lambda*/, const MaterialState* /*state*/, const PhotonPacket* /*pp*/) const
{
    return 0.;
}

////////////////////////////////////////////////////////////////////

double ElectronMix::opacitySca(double lambda, const MaterialState* state, const PhotonPacket* /*pp*/) const
{
    return state->numberDensity() * sectionSca(lambda);
}

////////////////////////////////////////////////////////////////////

double ElectronMix::opacityExt(double lambda, const MaterialState* state, const PhotonPacket* /*pp*/) const
{
    return state->numberDensity() * sectionSca(lambda);
}

////////////////////////////////////////////////////////////////////

void ElectronMix::peeloffScattering(double& I, double& Q, double& U, double& V, double& lambda, double w,
                                    Direction bfkobs, Direction bfky, const MaterialState* /*state*/,
                                    const PhotonPacket* pp) const
{
    if (_hasCompton && lambda < comptonWL)
        _cpf.peeloffScattering(I, lambda, w, pp->direction(), bfkobs);
    else
        _dpf.peeloffScattering(I, Q, U, V, w, pp->direction(), bfkobs, bfky, pp);
}

////////////////////////////////////////////////////////////////////

void ElectronMix::performScattering(double lambda, const MaterialState* state, PhotonPacket* pp) const
{
    // determine the new propagation direction, and if required,
    // update the wavelength or the polarization state of the photon packet
    Direction bfknew = (_hasCompton && lambda < comptonWL) ? _cpf.performScattering(lambda, pp->direction())
                                                           : _dpf.performScattering(pp->direction(), pp);

    // execute the scattering event in the photon packet
    pp->scatter(bfknew, state->bulkVelocity(), lambda);
}

////////////////////////////////////////////////////////////////////
