/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "DipolePhaseFunction.hpp"
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

void DipolePhaseFunction::initialize(Random* random, bool includePolarization)
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

double DipolePhaseFunction::phaseFunctionValueForCosine(double costheta) const
{
    return 0.75 * (costheta * costheta + 1.);
}

////////////////////////////////////////////////////////////////////

double DipolePhaseFunction::generateCosineFromPhaseFunction() const
{
    double X = _random->uniform();
    double p = cbrt(4. * X - 2. + sqrt(16. * X * (X - 1.) + 5.));
    return p - 1. / p;
}

////////////////////////////////////////////////////////////////////

double DipolePhaseFunction::phaseFunctionValue(double theta, double phi, const StokesVector* sv) const
{
    double costheta = cos(theta);
    double S11 = costheta * costheta + 1.;
    double S12 = costheta * costheta - 1.;
    double polDegree = sv->linearPolarizationDegree();
    double polAngle = sv->polarizationAngle();
    return 0.75 * (S11 + polDegree * S12 * cos(2. * (phi - polAngle)));
}

////////////////////////////////////////////////////////////////////

std::pair<double, double> DipolePhaseFunction::generateAnglesFromPhaseFunction(const StokesVector* sv) const
{
    // get the theta cosine just as in the unpolarized case
    double costheta = generateCosineFromPhaseFunction();

    // construct and sample from the normalized cumulative distribution of phi for this wavelength and theta angle
    double polDegree = sv->linearPolarizationDegree();
    double polAngle = sv->polarizationAngle();
    double PF = polDegree * (costheta * costheta + 1.) / (costheta * costheta - 1.) / (4 * M_PI);
    double cos2polAngle = cos(2 * polAngle) * PF;
    double sin2polAngle = sin(2 * polAngle) * PF;
    double phi = _random->cdfLinLin(_phiv, _phi1v + cos2polAngle * _phisv + sin2polAngle * _phicv);

    // return the result
    return std::make_pair(acos(costheta), phi);
}

////////////////////////////////////////////////////////////////////

void DipolePhaseFunction::applyMueller(double theta, StokesVector* sv) const
{
    double costheta = cos(theta);
    sv->applyMueller(costheta * costheta + 1., costheta * costheta - 1., 2. * costheta, 0.);
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

void DipolePhaseFunction::peeloffScattering(double& I, double& Q, double& U, double& V, Direction bfk, Direction bfkobs,
                                            Direction bfky, const StokesVector* sv) const
{
    if (!_includePolarization)
    {
        // calculate the value of the material-specific phase function
        double costheta = Vec::dot(bfk, bfkobs);
        double value = phaseFunctionValueForCosine(costheta);

        // accumulate the weighted sum in the intensity
        I += value;
    }
    else
    {
        // calculate the value of the material-specific phase function
        double theta = acos(Vec::dot(bfk, bfkobs));
        double phi = angleBetweenScatteringPlanes(sv->normal(), bfk, bfkobs);
        double value = phaseFunctionValue(theta, phi, sv);

        // copy the polarization state so we can change it without affecting the incoming stokes vector
        StokesVector svnew = *sv;

        // rotate the Stokes vector reference direction into the scattering plane
        svnew.rotateIntoPlane(bfk, bfkobs);

        // apply the Mueller matrix
        applyMueller(theta, &svnew);

        // rotate the Stokes vector reference direction parallel to the instrument frame y-axis
        // it is given bfkobs because the photon is at this point aimed towards the observer
        svnew.rotateIntoPlane(bfkobs, bfky);

        // acumulate the weighted sum of all Stokes components to support polarization
        I += value * svnew.stokesI();
        Q += value * svnew.stokesQ();
        U += value * svnew.stokesU();
        V += value * svnew.stokesV();
    }
}

////////////////////////////////////////////////////////////////////

Direction DipolePhaseFunction::performScattering(Direction bfk, StokesVector* sv) const
{
    // determine the new propagation direction, and if required, update the polarization state of the photon packet
    if (!_includePolarization)
    {
        // sample a scattering angle from the dipole phase function
        double costheta = generateCosineFromPhaseFunction();
        return _random->direction(bfk, costheta);
    }
    else
    {
        // sample the angles between the previous and new direction from the dipole phase function,
        // given the incoming polarization state
        double theta, phi;
        std::tie(theta, phi) = generateAnglesFromPhaseFunction(sv);

        // rotate the Stokes vector (and the scattering plane) of the photon packet
        sv->rotateStokes(phi, bfk);

        // apply the Mueller matrix to the Stokes vector of the photon packet
        applyMueller(theta, sv);

        // rotate the propagation direction in the scattering plane
        Vec newdir = bfk * cos(theta) + Vec::cross(sv->normal(), bfk) * sin(theta);

        // normalize the new direction to prevent degradation
        return Direction(newdir / newdir.norm());
    }
}

////////////////////////////////////////////////////////////////////
