/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "SingleGrainDustMix.hpp"
#include "NR.hpp"
#include "PhotonPacket.hpp"
#include "Random.hpp"

////////////////////////////////////////////////////////////////////

namespace
{
    // arbitrary constants determining the resolution for disretizing the scattering angles
    constexpr int Ntheta = 181;                 // theta from 0 to pi, index t
    constexpr double dtheta = M_PI/(Ntheta-1);
    constexpr int Nphi = 361;                   // phi from 0 to 2 pi, index f
    constexpr double dphi = 2*M_PI/(Nphi-1);
}

////////////////////////////////////////////////////////////////////

void SingleGrainDustMix::setupSelfBefore()
{
    MaterialMix::setupSelfBefore();

    // get the basic optical properties
    string name = resourceNameForOpticalProps();
    _sigmaabs.open(this, name, "lambda(m)", "sigmaabs(m2/H)");
    _sigmasca.open(this, name, "lambda(m)", "sigmasca(m2/H)");
    _asymmpar.open(this, name, "lambda(m)", "g(1)");
    _mu = StoredTable<1>(this, name, "lambda(m)", "mu(kg/H)")[1.];  // get mu value for arbitrary wavelength

    // get the Mueller matrix elements
    name = resourceNameForMuellerMatrix();
    if (!name.empty())
    {
        _S11.open(this, name, "lambda(m),theta(rad)", "S11(1)");
        _S12.open(this, name, "lambda(m),theta(rad)", "S12(1)");
        _S33.open(this, name, "lambda(m),theta(rad)", "S33(1)");
        _S34.open(this, name, "lambda(m),theta(rad)", "S34(1)");
        _polarization = true;
    }

    // create a table containing the theta value corresponding to each index in our discretization
    _thetav.resize(Ntheta);
    for (int t=0; t<Ntheta; t++) _thetav[t] = t * dtheta;

    // create tables listing phi, phi/(2 pi), sin(2 phi) and 1-cos(2 phi) for a number of phi indices
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

double SingleGrainDustMix::sigmaabs(double lambda) const
{
    return _sigmaabs[lambda];
}

////////////////////////////////////////////////////////////////////

double SingleGrainDustMix::sigmasca(double lambda) const
{
    return _sigmasca[lambda];
}

////////////////////////////////////////////////////////////////////

double SingleGrainDustMix::sigmaext(double lambda) const
{
    return _sigmaabs[lambda] + _sigmasca[lambda];
}

////////////////////////////////////////////////////////////////////

double SingleGrainDustMix::kappaabs(double lambda) const
{
    return _sigmaabs[lambda] / _mu;
}

////////////////////////////////////////////////////////////////////

double SingleGrainDustMix::kappasca(double lambda) const
{
    return _sigmasca[lambda] / _mu;
}

////////////////////////////////////////////////////////////////////

double SingleGrainDustMix::kappaext(double lambda) const
{
    return (_sigmaabs[lambda] + _sigmasca[lambda]) / _mu;
}

////////////////////////////////////////////////////////////////////

double SingleGrainDustMix::albedo(double lambda) const
{
    double sigmaabs = _sigmaabs[lambda];
    double sigmasca = _sigmasca[lambda];
    double sigmaext = sigmaabs + sigmasca;
    return sigmaext > 0. ? sigmasca/sigmaext : 0.;
}

////////////////////////////////////////////////////////////////////

bool SingleGrainDustMix::hasScatteringPolarization() const
{
    return _polarization;
}

////////////////////////////////////////////////////////////////////

Direction SingleGrainDustMix::scatteringDirectionAndPolarization(StokesVector* out, const PhotonPacket* pp) const
{
    if (_polarization)
    {
        // determine the angles between the previous and new direction
        double lambda = pp->wavelength();
        double theta = sampleTheta(lambda);
        double phi = samplePhi(lambda, theta, pp->linearPolarizationDegree(), pp->polarizationAngle());

        // rotate the Stokes vector (and the scattering plane)
        *out = *pp;
        out->rotateStokes(phi, pp->direction());

        // apply Mueller matrix
        out->applyMueller(_S11(lambda,theta), _S12(lambda,theta), _S33(lambda,theta), _S34(lambda,theta));

        // rotate the propagation direction in the scattering plane
        Vec newdir = pp->direction()*cos(theta) + Vec::cross(out->normal(), pp->direction())*sin(theta);

        // normalize direction to prevent degradation
        return Direction(newdir/newdir.norm());
    }
    else
    {
        double g = _asymmpar[pp->wavelength()];
        if (fabs(g)<1e-6) return random()->direction();
        double f = ((1.0-g)*(1.0+g))/(1.0-g+2.0*g*random()->uniform());
        double costheta = (1.0+g*g-f*f)/(2.0*g);
        return random()->direction(pp->direction(), costheta);
    }
}

////////////////////////////////////////////////////////////////////

void SingleGrainDustMix::scatteringPeelOffPolarization(StokesVector* out, const PhotonPacket* pp, Direction bfknew,
                                                       Direction /*bfkx*/, Direction bfky)
{
    if (_polarization)
    {
        // copy the polarization state
        *out = *pp;

        // rotate the Stokes vector reference direction into the scattering plane
        out->rotateIntoPlane(pp->direction(),bfknew);

        // apply the Mueller matrix
        double lambda = pp->wavelength();
        double theta = acos(Vec::dot(pp->direction(),bfknew));
        out->applyMueller(_S11(lambda,theta), _S12(lambda,theta), _S33(lambda,theta), _S34(lambda,theta));

        // rotate the Stokes vector reference direction parallel to the instrument frame y-axis
        // it is given bfknew, because the photon is at this point aimed towards the observer,
        // but the propagation direction has not been updated.
        out->rotateIntoPlane(bfknew,bfky);
    }
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

double SingleGrainDustMix::phaseFunctionValue(const PhotonPacket* pp, Direction bfknew) const
{
    if (_polarization)
    {
        // calculate the phase function normalization factor for the photon packet wavelength
        double lambda = pp->wavelength();
        double sum = 0.;
        for (int t=1; t<Ntheta-1; t++)        // sin(0) and sin(pi) are zero anyway
        {
            double theta = _thetav[t];
            sum += _S11(lambda,theta)*sin(theta)*dtheta;
        }
        double normalization = 2.0/sum;

        // determine the scattering angles
        double phi = angleBetweenScatteringPlanes(pp->normal(), pp->direction(), bfknew);
        double theta = acos(Vec::dot(pp->direction(),bfknew));

        // calculate the phase function value
        double polDegree = pp->linearPolarizationDegree();
        double polAngle = pp->polarizationAngle();
        return normalization * (_S11(lambda,theta) + polDegree*_S12(lambda,theta)*cos(2.*(phi-polAngle)));
    }
    else
    {
        double cosalpha = Direction::dot(pp->direction(), bfknew);
        double g = _asymmpar[pp->wavelength()];
        double t = 1.0+g*g-2*g*cosalpha;
        return (1.0-g)*(1.0+g)/sqrt(t*t*t);
    }
}

////////////////////////////////////////////////////////////////////

double SingleGrainDustMix::sampleTheta(double lambda) const
{
    // construct the normalized cumulative distribution of theta for this wavelength
    Array thetaXv;
    NR::cdf(thetaXv, Ntheta-1, [this,lambda](int t){ return _S11(lambda,_thetav[t+1])*sin(_thetav[t+1])*dtheta; });

    // sample theta from it
    return random()->cdfLinLin(_thetav, thetaXv);
}

////////////////////////////////////////////////////////////////////

double SingleGrainDustMix::samplePhi(double lambda, double theta, double polDegree, double polAngle) const
{
    // construct and sample from the normalized cumulative distribution of phi for this wavelength and theta angle
    double PF = polDegree * _S12(lambda,theta)/_S11(lambda,theta) / (4*M_PI);
    double cos2polAngle = cos(2*polAngle) * PF;
    double sin2polAngle = sin(2*polAngle) * PF;
    return random()->cdfLinLin(_phiv, _phi1v + cos2polAngle*_phisv + sin2polAngle*_phicv);
}

////////////////////////////////////////////////////////////////////
