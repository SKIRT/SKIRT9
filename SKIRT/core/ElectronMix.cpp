/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "ElectronMix.hpp"
#include "Constants.hpp"
#include "PhotonPacket.hpp"

////////////////////////////////////////////////////////////////////

double ElectronMix::sigmaabs(double /*lambda*/) const
{
    return 0.;
}

////////////////////////////////////////////////////////////////////

double ElectronMix::sigmasca(double /*lambda*/) const
{
    return Constants::sigmaThomson();
}

////////////////////////////////////////////////////////////////////

double ElectronMix::sigmaext(double /*lambda*/) const
{
    return Constants::sigmaThomson();
}

////////////////////////////////////////////////////////////////////

double ElectronMix::albedo(double /*lambda*/) const
{
    return 1.;
}

////////////////////////////////////////////////////////////////////

bool ElectronMix::hasScatteringPolarization() const
{
    return true;
}

////////////////////////////////////////////////////////////////////

Direction ElectronMix::scatteringDirectionAndPolarization(StokesVector* out, const PhotonPacket* pp) const
{
    // determine the angles between the previous and new direction
    double theta = sampleTheta();
    double costheta = cos(theta);
    double sintheta = sin(theta);
    double phi = samplePhi(costheta, pp->linearPolarizationDegree(), pp->polarizationAngle());

    // rotate the Stokes vector (and the scattering plane)
    *out = *pp;
    out->rotateStokes(phi, pp->direction());

    // apply Mueller matrix
    out->applyMueller(costheta*costheta + 1., costheta*costheta - 1., 2.*costheta, 0.);

    // rotate the propagation direction in the scattering plane
    Vec newdir = pp->direction()*costheta + Vec::cross(out->normal(), pp->direction())*sintheta;

    // normalize direction to prevent degradation
    return Direction(newdir/newdir.norm());
}

////////////////////////////////////////////////////////////////////

void ElectronMix::scatteringPeelOffPolarization(StokesVector* out, const PhotonPacket* pp, Direction bfknew,
                                                Direction /*bfkx*/, Direction bfky)
{
    // copy the polarization state
    *out = *pp;

    // rotate the Stokes vector reference direction into the scattering plane
    out->rotateIntoPlane(pp->direction(),bfknew);

    // apply the Mueller matrix
    double costheta = Vec::dot(pp->direction(),bfknew);
    out->applyMueller(costheta*costheta + 1., costheta*costheta - 1., 2.*costheta, 0.);

    // rotate the Stokes vector reference direction parallel to the instrument frame y-axis
    // it is given bfknew, because the photon is at this point aimed towards the observer,
    // but the propagation direction has not been updated.
    out->rotateIntoPlane(bfknew,bfky);
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
        Vec nc = Vec::cross(kc,kn);
        nc /= nc.norm();
        double cosphi = Vec::dot(np,nc);
        double sinphi = Vec::dot(Vec::cross(np,nc), kc);
        double phi = atan2(sinphi,cosphi);
        if (std::isfinite(phi)) return phi;
        return 0;
    }
}

////////////////////////////////////////////////////////////////////

double ElectronMix::phaseFunctionValue(const PhotonPacket* pp, Direction bfknew) const
{
    // determine the scattering angles
    double phi = angleBetweenScatteringPlanes(pp->normal(), pp->direction(), bfknew);
    double costheta = Vec::dot(pp->direction(),bfknew);

    // calculate the phase function value
    double normalization = 0.75;
    double S11 = costheta*costheta + 1.;
    double S12 = costheta*costheta - 1.;
    double polDegree = pp->linearPolarizationDegree();
    double polAngle = pp->polarizationAngle();
    return normalization * (S11 + polDegree*S12*cos(2.*(phi-polAngle)));
}

////////////////////////////////////////////////////////////////////
