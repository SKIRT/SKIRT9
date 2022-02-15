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
    const double comptonWL = 9.999999e-9;  // 10 nm
}

////////////////////////////////////////////////////////////////////

void ElectronMix::setupSelfBefore()
{
    MaterialMix::setupSelfBefore();

    // determine whether velocity dispersion is enabled
    _hasDispersion = includeThermalDispersion() && !find<Configuration>()->oligochromatic();

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

bool ElectronMix::hasScatteringDispersion() const
{
    // this function may be called before setup() has been invoked, so we need to figure out the result ourselves

    // determine whether velocity dispersion is enabled
    bool hasDispersion = includeThermalDispersion() && !find<Configuration>()->oligochromatic();

    // determine whether we need Compton scattering because the simulation may have short wavelengths
    if (!hasDispersion && !includePolarization())
    {
        auto range = find<Configuration>()->simulationWavelengthRange();
        hasDispersion = range.min() < comptonWL;
    }

    return hasDispersion;
}

////////////////////////////////////////////////////////////////////

vector<StateVariable> ElectronMix::specificStateVariableInfo() const
{
    vector<StateVariable> result{StateVariable::numberDensity()};
    if (_hasDispersion) result.push_back(StateVariable::temperature());
    return result;
}

////////////////////////////////////////////////////////////////////

void ElectronMix::initializeSpecificState(MaterialState* state, double /*metallicity*/, double temperature,
                                          const Array& /*params*/) const
{
    // leave the temperature at zero if the cell does not contain any material for this component
    if (_hasDispersion && state->numberDensity() > 0.)
    {
        // if no temperature was imported, use default value
        if (temperature < 0) temperature = defaultTemperature();

        // make sure the temperature is at least the local universe CMB temperature
        state->setTemperature(max(Constants::Tcmb(), temperature));
    }
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

namespace
{
    // draw a random electron velocity unless a previous peel-off stored this already
    void generateElectronVelocityIfNeeded(PhotonPacket::ScatteringInfo* scatinfo, double T, Random* random)
    {
        if (!scatinfo->valid)
        {
            scatinfo->valid = true;

            // for high temperatures the generated velocity can be relativistic or even above the speed of light,
            // in which case our non-relativistic Doppler shift formulas produce negative wavelengths;
            // we thus reject any velocities above c/3
            constexpr double vmax2 = Constants::c() * Constants::c() / 9.;
            while (true)
            {
                scatinfo->velocity = sqrt(Constants::k() / Constants::Melectron() * T) * random->maxwell();
                if (scatinfo->velocity.norm2() < vmax2) break;
            }
        }
    }
}

////////////////////////////////////////////////////////////////////

void ElectronMix::peeloffScattering(double& I, double& Q, double& U, double& V, double& lambda, Direction bfkobs,
                                    Direction bfky, const MaterialState* state, const PhotonPacket* pp) const
{
    // if we have dispersion, adjust the incoming wavelength to the electron rest frame
    PhotonPacket::ScatteringInfo* scatinfo = nullptr;
    if (_hasDispersion)
    {
        // draw a random electron velocity unless a previous peel-off stored this already
        scatinfo = const_cast<PhotonPacket*>(pp)->getScatteringInfo();
        generateElectronVelocityIfNeeded(scatinfo, state->temperature(), random());

        // adjust the wavelength
        lambda = PhotonPacket::shiftedReceptionWavelength(lambda, pp->direction(), scatinfo->velocity);
    }

    // perform the scattering event in the electron rest frame
    if (_hasCompton && lambda < comptonWL)
        _cpf.peeloffScattering(I, lambda, pp->direction(), bfkobs);
    else
        _dpf.peeloffScattering(I, Q, U, V, pp->direction(), bfkobs, bfky, pp);

    // if we have dispersion, adjust the outgoing wavelength from the electron rest frame
    if (_hasDispersion)
    {
        lambda = PhotonPacket::shiftedEmissionWavelength(lambda, bfkobs, scatinfo->velocity);
    }
}

////////////////////////////////////////////////////////////////////

void ElectronMix::performScattering(double lambda, const MaterialState* state, PhotonPacket* pp) const
{
    // if we have dispersion, adjust the incoming wavelength to the electron rest frame
    PhotonPacket::ScatteringInfo* scatinfo = nullptr;
    if (_hasDispersion)
    {
        // draw a random electron velocity unless a previous peel-off stored this already
        scatinfo = pp->getScatteringInfo();
        generateElectronVelocityIfNeeded(scatinfo, state->temperature(), random());

        // adjust the wavelength
        lambda = PhotonPacket::shiftedReceptionWavelength(lambda, pp->direction(), scatinfo->velocity);
    }

    // determine the new propagation direction, and if required,
    // update the wavelength or the polarization state of the photon packet
    Direction bfknew = (_hasCompton && lambda < comptonWL) ? _cpf.performScattering(lambda, pp->direction())
                                                           : _dpf.performScattering(pp->direction(), pp);

    // if we have dispersion, adjust the outgoing wavelength from the electron rest frame
    if (_hasDispersion)
    {
        lambda = PhotonPacket::shiftedEmissionWavelength(lambda, bfknew, scatinfo->velocity);
    }

    // execute the scattering event in the photon packet
    pp->scatter(bfknew, state->bulkVelocity(), lambda);
}

////////////////////////////////////////////////////////////////////

double ElectronMix::indicativeTemperature(const MaterialState* state, const Array& /*Jv*/) const
{
    return _hasDispersion ? state->temperature() : 0.;
}

////////////////////////////////////////////////////////////////////
