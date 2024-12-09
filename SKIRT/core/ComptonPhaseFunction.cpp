/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "ComptonPhaseFunction.hpp"
#include "Constants.hpp"
#include "NR.hpp"
#include "Random.hpp"
#include "StokesVector.hpp"

////////////////////////////////////////////////////////////////////

namespace
{
    // arbitrary constants determining the resolution for discretizing the scattering angles
    constexpr int Nphi = 361;  // phi from 0 to 2 pi, index f
    constexpr double dphi = 2 * M_PI / (Nphi - 1);
}

////////////////////////////////////////////////////////////////////

namespace
{
    // returns the photon energy scaled to the electron rest energy: h nu / m_e c^2
    double scaledEnergy(double lambda)
    {
        return (Constants::h() / Constants::Melectron() / Constants::c()) / lambda;
    }

    // returns the Compton scattering cross section (relative to the Thomson cross section) for a given scaled energy
    double comptonSection(double x)
    {
        double x2 = x * x;
        double x3 = x2 * x;
        double x1 = 1. + x;
        double x12 = 1. + 2. * x;
        return 0.75 * (x1 / (x12 * x12) + 2. / x2 + (1. / (2. * x) - x1 / x3) * log(x12));
    }

    // returns the inverse Compton factor for a given scaled energy and scattering angle cosine
    double inverseComptonfactor(double x, double costheta)
    {
        return 1 + x * (1 - costheta);
    }

    // returns the Compton factor for a given scaled energy and scattering angle cosine
    double comptonFactor(double x, double costheta)
    {
        return 1. / inverseComptonfactor(x, costheta);
    }
}

////////////////////////////////////////////////////////////////////

void ComptonPhaseFunction::initialize(Random* random, bool includePolarization)
{
    // cache random number generator and polarization flag
    _random = random;
    _includePolarization = includePolarization;

    if (includePolarization)
    {
        // create tables listing phi, phi/(2 pi), sin(2 phi) and 1-cos(2 phi) for a number of f indices
        _phiv.resize(Nphi);
        _phi1v.resize(Nphi);
        _phisv.resize(Nphi);
        _phicv.resize(Nphi);
        for (int f = 0; f < Nphi; f++)
        {
            double phi = f * dphi;
            _phiv[f] = phi;
            _phi1v[f] = phi / (2 * M_PI);
            _phisv[f] = sin(2 * phi);
            _phicv[f] = 1 - cos(2 * phi);
        }
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
    double sin2theta = (1. - costheta) * (1. + costheta);
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

double ComptonPhaseFunction::phaseFunctionValue(double x, double costheta, double phi, const StokesVector* sv) const
{
    double C = comptonFactor(x, costheta);
    double sin2theta = (1. - costheta) * (1. + costheta);
    double S12 = -C * C * sin2theta;
    double S11 = C * C * C + C + S12;
    double polDegree = sv->linearPolarizationDegree();
    double polAngle = sv->polarizationAngle();
    return 0.75 / comptonSection(x) * (S11 + polDegree * S12 * cos(2. * (phi - polAngle)));
}

////////////////////////////////////////////////////////////////////

double ComptonPhaseFunction::generateAzimuthFromPhaseFunction(double x, const StokesVector* sv, double costheta) const
{
    // construct and sample from the normalized cumulative distribution of phi for this wavelength and theta angle
    double C = comptonFactor(x, costheta);
    double sin2theta = (1. - costheta) * (1. + costheta);
    double S12 = -C * C * sin2theta;
    double S11 = C * C * C + C + S12;
    double polDegree = sv->linearPolarizationDegree();
    double polAngle = sv->polarizationAngle();
    double PF = polDegree * S12 / S11 / (4 * M_PI);
    double cos2polAngle = cos(2 * polAngle) * PF;
    double sin2polAngle = sin(2 * polAngle) * PF;
    double phi = _random->cdfLinLin(_phiv, _phi1v + cos2polAngle * _phisv + sin2polAngle * _phicv);
    return phi;
}

////////////////////////////////////////////////////////////////////

void ComptonPhaseFunction::applyMueller(double x, double costheta, StokesVector* sv) const
{
    double C = comptonFactor(x, costheta);
    double C2 = C * C;
    double C3 = C2 * C;
    double sin2theta = (1. - costheta) * (1. + costheta);
    double S12 = -C2 * sin2theta;
    double S11 = C3 + C + S12;
    double S22 = C2 * (1. + costheta * costheta);
    double S33 = 2. * C2 * costheta;
    double S44 = (C3 + C) * costheta;
    sv->applyMueller(S11, S12, S22, S33, 0., S44);
}

////////////////////////////////////////////////////////////////////

namespace
{
    // This helper function returns the angle phi between the previous and current scattering planes
    // given the normal to the previous scattering plane and the current and new propagation directions
    // of the photon packet. The function returns a zero angle if the light is unpolarized or when the
    // current scattering event is completely forward or backward.
    double angleBetweenScatteringPlanes(Direction np, Direction kc, Direction kn)
    {
        Vec nc = Vec::cross(kc, kn);
        nc /= nc.norm();
        double cosphi = Vec::dot(np, nc);
        double sinphi = Vec::dot(Vec::cross(np, nc), kc);
        double phi = atan2(sinphi, cosphi);
        if (std::isfinite(phi)) return phi;
        return 0.;
    }
}

////////////////////////////////////////////////////////////////////

void ComptonPhaseFunction::peeloffScattering(double& I, double& Q, double& U, double& V, double& lambda, Direction bfk,
                                             Direction bfkobs, Direction bfky, const StokesVector* sv) const
{
    // get the scaled energy and the scattering angle cosine
    double x = scaledEnergy(lambda);
    double costheta = Vec::dot(bfk, bfkobs);

    if (!_includePolarization)
    {
        // calculate the value of the phase function
        double value = phaseFunctionValueForCosine(x, costheta);

        // store this value as the intensity
        I = value;
    }
    else
    {
        // calculate the value of the material-specific phase function
        double phi = angleBetweenScatteringPlanes(sv->normal(), bfk, bfkobs);
        double value = phaseFunctionValue(x, costheta, phi, sv);

        // copy the polarization state so we can change it without affecting the incoming stokes vector
        StokesVector svnew = *sv;

        // rotate the Stokes vector reference direction into the scattering plane
        svnew.rotateIntoPlane(bfk, bfkobs);

        // apply the Mueller matrix
        applyMueller(x, costheta, &svnew);

        // rotate the Stokes vector reference direction parallel to the instrument frame y-axis
        // it is given bfkobs because the photon is at this point aimed towards the observer
        svnew.rotateIntoPlane(bfkobs, bfky);

        // store the new Stokes vector components
        I = value * svnew.stokesI();
        Q = value * svnew.stokesQ();
        U = value * svnew.stokesU();
        V = value * svnew.stokesV();
    }

    // adjust the wavelength
    lambda *= inverseComptonfactor(x, costheta);
}

////////////////////////////////////////////////////////////////////

Direction ComptonPhaseFunction::performScattering(double& lambda, Direction bfk, StokesVector* sv) const
{
    // get the scaled energy
    double x = scaledEnergy(lambda);
    double costheta = generateCosineFromPhaseFunction(x);

    // adjust the wavelength
    lambda *= inverseComptonfactor(x, costheta);

    // determine the new propagation direction, and if required, update the polarization state of the photon packet
    if (!_includePolarization)
    {
        return _random->direction(bfk, costheta);
    }
    else
    {
        // sample the azimuthal scattering angle, given the incoming polarization state amd the inclination cosine
        double phi = generateAzimuthFromPhaseFunction(x, sv, costheta);

        // rotate the Stokes vector (and the scattering plane) of the photon packet
        sv->rotateStokes(phi, bfk);

        // apply the Mueller matrix to the Stokes vector of the photon packet
        applyMueller(x, costheta, sv);

        // rotate the propagation direction in the scattering plane
        // (re)normalize the new direction to prevent degradation
        return Direction(bfk * costheta + Vec::cross(sv->normal(), bfk) * sin(acos(costheta)), true);
    }
}

////////////////////////////////////////////////////////////////////
