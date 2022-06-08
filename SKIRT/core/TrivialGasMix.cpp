/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "TrivialGasMix.hpp"
#include "Constants.hpp"
#include "MaterialState.hpp"
#include "PhotonPacket.hpp"
#include "Random.hpp"

////////////////////////////////////////////////////////////////////

MaterialMix::MaterialType TrivialGasMix::materialType() const
{
    return MaterialType::Gas;
}

////////////////////////////////////////////////////////////////////

bool TrivialGasMix::hasNegativeExtinction() const
{
    // capture the border case where the magnitudes are nonzero and equal
    return absorptionCrossSection() < 0. && (absorptionCrossSection() + scatteringCrossSection()) <= 0.;
}

////////////////////////////////////////////////////////////////////

vector<StateVariable> TrivialGasMix::specificStateVariableInfo() const
{
    return vector<StateVariable>{StateVariable::numberDensity()};
}

////////////////////////////////////////////////////////////////////

double TrivialGasMix::mass() const
{
    return Constants::Mproton();
}

////////////////////////////////////////////////////////////////////

double TrivialGasMix::sectionAbs(double /*lambda*/) const
{
    return absorptionCrossSection();
}

////////////////////////////////////////////////////////////////////

double TrivialGasMix::sectionSca(double /*lambda*/) const
{
    return scatteringCrossSection();
}

////////////////////////////////////////////////////////////////////

double TrivialGasMix::sectionExt(double /*lambda*/) const
{
    return absorptionCrossSection() + scatteringCrossSection();
}

////////////////////////////////////////////////////////////////////

double TrivialGasMix::asymmpar(double /*lambda*/) const
{
    return asymmetryParameter();
}

////////////////////////////////////////////////////////////////////

double TrivialGasMix::opacityAbs(double /*lambda*/, const MaterialState* state, const PhotonPacket* /*pp*/) const
{
    return absorptionCrossSection() * state->numberDensity();
}

////////////////////////////////////////////////////////////////////

double TrivialGasMix::opacitySca(double /*lambda*/, const MaterialState* state, const PhotonPacket* /*pp*/) const
{
    return scatteringCrossSection() * state->numberDensity();
}

////////////////////////////////////////////////////////////////////

double TrivialGasMix::opacityExt(double /*lambda*/, const MaterialState* state, const PhotonPacket* /*pp*/) const
{
    return (absorptionCrossSection() + scatteringCrossSection()) * state->numberDensity();
}

////////////////////////////////////////////////////////////////////

void TrivialGasMix::peeloffScattering(double& I, double& /*Q*/, double& /*U*/, double& /*V*/, double& /*lambda*/,
                                      Direction bfkobs, Direction /*bfky*/, const MaterialState* /*state*/,
                                      const PhotonPacket* pp) const
{
    // calculate the value of the Henyey-Greenstein phase function
    double costheta = Vec::dot(pp->direction(), bfkobs);
    double g = asymmetryParameter();
    double t = 1. + g * g - 2. * g * costheta;
    double value = (1. - g) * (1. + g) / sqrt(t * t * t);

    // accumulate the weighted sum in the intensity
    I += value;
}

////////////////////////////////////////////////////////////////////

void TrivialGasMix::performScattering(double lambda, const MaterialState* state, PhotonPacket* pp) const
{
    // determine the new propagation direction
    // sample a scattering angle from the Henyey-Greenstein phase function
    // handle isotropic scattering separately because the HG sampling procedure breaks down in this case
    Direction bfknew;
    double g = asymmetryParameter();
    if (fabs(g) < 1e-6)
    {
        bfknew = random()->direction();
    }
    else
    {
        double f = ((1.0 - g) * (1.0 + g)) / (1.0 - g + 2.0 * g * random()->uniform());
        double costheta = (1.0 + g * g - f * f) / (2.0 * g);
        bfknew = random()->direction(pp->direction(), costheta);
    }

    // execute the scattering event in the photon packet
    pp->scatter(bfknew, state->bulkVelocity(), lambda);
}

////////////////////////////////////////////////////////////////////

double TrivialGasMix::indicativeTemperature(const MaterialState* /*state*/, const Array& /*Jv*/) const
{
    return 0.;
}

////////////////////////////////////////////////////////////////////
