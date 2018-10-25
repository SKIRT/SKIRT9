/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "ElectronMix.hpp"
#include "Constants.hpp"
#include "StokesVector.hpp"
#include "Random.hpp"

////////////////////////////////////////////////////////////////////

namespace
{
    // arbitrary constants determining the resolution for disretizing the scattering angles
    constexpr int Nphi = 361;                   // phi from 0 to 2 pi, index f
    constexpr double dphi = 2*M_PI/(Nphi-1);
}

////////////////////////////////////////////////////////////////////

void ElectronMix::setupSelfBefore()
{
    MaterialMix::setupSelfBefore();

    // create tables listing phi, phi/(2 pi), sin(2 phi) and 1-cos(2 phi) for a number of f indices
    _phiv.resize(Nphi);
    _phi1v.resize(Nphi);
    _phisv.resize(Nphi);
    _phicv.resize(Nphi);
    for (int f=0; f<Nphi; f++)
    {
        double phi = f * dphi;
        _phiv[f] = phi;
        _phi1v[f] = phi/(2*M_PI);
        _phisv[f] = sin(2*phi);
        _phicv[f] = 1-cos(2*phi);
    }
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

double ElectronMix::sectionExtSelf(double /*lambda*/) const
{
    return Constants::sigmaThomson();
}

////////////////////////////////////////////////////////////////////

double ElectronMix::albedo(double /*lambda*/) const
{
    return 1.;
}

////////////////////////////////////////////////////////////////////

double ElectronMix::phaseFunctionValueForCosine(double /*lambda*/, double costheta) const
{
    return 0.75 * (costheta*costheta + 1.);
}

////////////////////////////////////////////////////////////////////

double ElectronMix::generateCosineFromPhaseFunction(double /*lambda*/) const
{
    double X = random()->uniform();
    double p = cbrt( 4.*X - 2. + sqrt(16.*X*(X-1.) + 5.) );
    return p - 1./p;
}

////////////////////////////////////////////////////////////////////

double ElectronMix::phaseFunctionValue(double /*lambda*/, double theta, double phi, const StokesVector* sv) const
{
    double costheta = cos(theta);
    double S11 = costheta*costheta + 1.;
    double S12 = costheta*costheta - 1.;
    double polDegree = sv->linearPolarizationDegree();
    double polAngle = sv->polarizationAngle();
    return 0.75 * (S11 + polDegree*S12*cos(2.*(phi-polAngle)));
}

////////////////////////////////////////////////////////////////////

std::pair<double, double> ElectronMix::generateAnglesFromPhaseFunction(double lambda, const StokesVector* sv) const
{
    // get the theta cosine just as in the unpolarized case
    double costheta = generateCosineFromPhaseFunction(lambda);

    // construct and sample from the normalized cumulative distribution of phi for this wavelength and theta angle
    double polDegree = sv->linearPolarizationDegree();
    double polAngle = sv->polarizationAngle();
    double PF = polDegree * (costheta*costheta + 1.) / (costheta*costheta - 1.) / (4*M_PI);
    double cos2polAngle = cos(2*polAngle) * PF;
    double sin2polAngle = sin(2*polAngle) * PF;
    double phi = random()->cdfLinLin(_phiv, _phi1v + cos2polAngle*_phisv + sin2polAngle*_phicv);

    // return the result
    return std::make_pair(acos(costheta), phi);
}

////////////////////////////////////////////////////////////////////

void ElectronMix::applyMueller(double /*lambda*/, double theta, StokesVector* sv) const
{
    double costheta = cos(theta);
    sv->applyMueller(costheta*costheta + 1., costheta*costheta - 1., 2.*costheta, 0.);
}

////////////////////////////////////////////////////////////////////
